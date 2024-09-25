nextflow.enable.dsl=2

params.clinvar="test_fixtures/chr18_clinvar.vcf"

params.gencode="test_fixtures/chr18_gencode.gff"

process vcf2bed {
  
  penv 'smp'
  cpus 1
  conda 'bedops=2.4.41 gawk'
  input:
    file clinvar_vcf

  output:
    file 'clinvar_bed.bed'

  shell:
  '''
   awk 'BEGIN { OFS="\t" } /^#/ { next; } { if($1 !~ /^chr/) { $1 = "chr"$1; } if($1 == "chrMT") { $1 = "chrM"; } print $1, $2-1, $2, "id="$3";ref="$4";alt="$5";qual="$6";filter="$7";"$8 }' !{clinvar_vcf} > clinvar_bed.bed
  '''
}

process intersectClinvarGencode {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 gawk'

        input:
          file clinvar_bed
          file gencode_gff

        output:
          file 'clinvar_regions.gff'
        shell:
        '''
        bedtools intersect -u -a !{gencode_gff} -b !{clinvar_bed} | awk '$3 == "gene" { print }' > clinvar_regions.gff
        '''
}


process gff2bed {

        penv 'smp'
        cpus 1
        conda 'bedops=2.4.41'

        input:
          file clinvar_gff

        output:
          file 'clinvar_regions.bed'
        shell:
        '''
        gff2bed < !{clinvar_gff} > clinvar_regions.bed
        '''
}

process sortBed {

        penv 'smp'
        cpus 1
        conda 'coreutils'
        publishDir 'clinvar_regions', mode: 'move'

        input:
          file clinvar_bed

        output:
          file 'clinvar_regions.sorted.bed'
        shell:
        '''
        sort -k1,1 -k2,2n !{clinvar_bed} > clinvar_regions.sorted.bed
        '''
}

workflow {
  clinvar_vcf = file(params.clinvar)

  gencode_gff = file(params.gencode)

  vcf2bed(clinvar_vcf)
  intersectClinvarGencode(vcf2bed.out, gencode_gff)
  gff2bed(intersectClinvarGencode.out)
  sortBed(gff2bed.out)
}

