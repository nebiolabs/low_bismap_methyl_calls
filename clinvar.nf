nextflow.enable.dsl=2

params.clinvar="test_fixtures/chr18_clinvar.vcf"

params.gencode="test_fixtures/chr18_gencode.gff"

process vcf2bed {
  
  input:
    file clinvar_vcf

  output:
    file 'clinvar_bed.bed'

  shell:
  '''
  vcf2bed < !{clinvar_vcf} | sed "s/^MT/M/g" | sed -E "s/(^[0-9]+)/chr\\1/g" > clinvar_bed.bed
  '''
}


process gff2bed {

  input:
    file gencode_gff

  output:
    file 'gencode_bed.bed'

  shell:
  '''
  gff2bed < !{gencode_gff} > gencode_bed.bed
  '''
}

process intersectClinvarGencode {

        publishDir 'clinvar_regions', mode: 'move'

        input:
          file clinvar_bed
          file gencode_bed

        output:
          file 'clinvar_regions.bed'
        shell:
        '''
        bedtools intersect -u -a !{gencode_bed} -b !{clinvar_bed} > clinvar_regions.bed
        '''
}

workflow {
  clinvar_vcf = file(params.clinvar)

  gencode_gff = file(params.gencode)

  vcf2bed(clinvar_vcf)
  gff2bed(gencode_gff)
  intersectClinvarGencode(vcf2bed.out, gff2bed.out)
}