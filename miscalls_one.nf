nextflow.enable.dsl=2

params.clinvarbed="test_fixtures/chr18_clinvar_regions.bed"

params.bismap="test_fixtures/chr18_bismap.bw"

params.bismap_bbm="test_fixtures/chr18_bismap.bbm"

params.ref_fai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.min_mapq=10

params.no_filter_mappability=false

params.bismap_cutoff=0.01

params.suffix="data"

process bigWigToBedGraph {

        conda 'ucsc-bigwigtobedgraph'

        input:
                file bismap_bw

        output:
                file 'bismap.bedGraph'

        shell:
        '''
        bigWigToBedGraph !{bismap_bw} bismap.bedGraph
        '''
}

process convertFaiToGenome {

        input:
                file ref_fai

        output:
                file 'ref.genome'

        shell:
        '''
        cut -f 1,2 !{ref_fai} > ref.genome
        '''
}

process extractBam {

        penv 'smp'
	cpus 8

        input:
          file ref
          val min_mapq
          each bam
          val filter_opt
          file ref_fai
          file bismap_bbm
	
        output:
          file '*_CpG.bedGraph'
		      file '*_CHG.bedGraph'
		      file '*_CHH.bedGraph'
		val bam
        shell:
        """
        /mnt/flash_scratch/ckumar/transfer/MethylDackel/MethylDackel extract -@ !{task.cpus} -q !{min_mapq} --CHH --CHG --OT 2,0,0,97 !{filter_opt} -o methyl_!{bam.split("/")[-1].replace(".bam\\n", "")} !{ref} !{bam} 
        """
}

process stringParsing {

        input:
		      val bam_raw
        
        output:
          stdout

   shell:
        """
        echo !{bam_raw.split("/")[-1].replace(".bam\n", "")}
        """
}



process combineCalls {

        conda 'coreutils'

        penv 'smp'
	cpus 8

        input:
                file methyl_cpg
		file methyl_chg
		file methyl_chh
		val n

        output:
                file 'sorted_noheader_*.bedGraph'
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{methyl_cpg} !{methyl_chg} !{methyl_chh} | head -n -3 > sorted_noheader_!{n.strip()}.bedGraph
        '''
}

process complementBedgraph {

        conda 'bedtools'

        input:
                file ref_genome
		file bismap_bedGraph_for_complement

        output:
                file 'bismap_zeroes.bedGraph'
        shell:
        '''
        bedtools complement -i !{bismap_bedGraph_for_complement} -g !{ref_genome} | awk '{ print $0"\t"0 }'> bismap_zeroes.bedGraph
        '''
}

process combineAllRegions {

        penv 'smp'
	cpus 8
        conda 'coreutils'

        input:

		file bismap_bedGraph_for_combine_all
		file zeroes_bedgraph

        output:
                file 'bismap_all.bedGraph'
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{bismap_bedGraph_for_combine_all} !{zeroes_bedgraph} > bismap_all.bedGraph
        '''
}

process intersectMin {

        conda 'bedtools'

        input:
                file clinvar_regions_for_min
		file bismap_bedgraph_for_min

        output:
                file 'clinvar_regions_min.bedGraph'
        shell:
        '''
        bedtools map -o min -c 4 -a !{clinvar_regions_for_min} -b !{bismap_bedgraph_for_min} > clinvar_regions_min.bedGraph
        '''
}

process intersectMean {

        conda 'bedtools'

        input:
                file clinvar_regions_for_mean
		file bismap_bedgraph_for_mean

        output:
                file 'clinvar_regions_mean.bedGraph'
        shell:
        '''
        bedtools map -o mean -c 4 -a !{clinvar_regions_for_mean} -b !{bismap_bedgraph_for_mean} > clinvar_regions_mean.bedGraph
        '''
}

process pasteMinMean {

        input:
                file clinvar_regions_mean_map
		file clinvar_regions_min_map

        output:
                file 'clinvar_regions_min_mean.bedGraph'
        shell:
        '''
        paste !{clinvar_regions_min_map} <(cut -f 11 !{clinvar_regions_mean_map}) > clinvar_regions_min_mean.bedGraph
        '''
}

process filterMinMean {

        input:
                file clinvar_regions_min_mean_map
                val bismap_cutoff

        output:
                file 'clinvar_regions_min_mean_low.bedGraph'
        shell:
        '''
        awk '$11 <= !{bismap_cutoff} { print }' !{clinvar_regions_min_mean_map} > clinvar_regions_min_mean_low.bedGraph
        '''
}

process gsortMinMean {

	cpus 8
        penv 'smp'
        conda 'coreutils'

        input:
                file clinvar_regions_min_mean_map_low_unsorted

        output:
                file 'clinvar_annotations_low.bedGraph'
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{clinvar_regions_min_mean_map_low_unsorted} > clinvar_annotations_low.bedGraph
        '''
}

process getLowRegions {

        input:

		file bismap_bedGraph_for_awk_low_regions

        output:
                file 'bismap_low.bedGraph'
        shell:
        '''
        awk '$4 <= 0.01 { print }' !{bismap_bedGraph_for_awk_low_regions} > bismap_low.bedGraph
        '''
}



process combineLowRegions {

	cpus 8
        penv 'smp'
        conda 'coreutils'

        input:

		file low_bedgraph
		file zeroes_bedgraph

        output:
                file 'bismap_low.bedGraph'
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{low_bedgraph} !{zeroes_bedgraph} > bismap_low.bedGraph
        '''
}

process intersectClinvarBismapLow {

        conda 'bedtools'

        input:
                file clinvar_annotations_low
		file bismap_low

        output:
		file 'clinvar_low_regions.bed'
        shell:
        '''
        bedtools intersect -sorted -a !{bismap_low} -b !{clinvar_annotations_low} > clinvar_low_regions.bed
        '''
}

process intersectCalls {

	publishDir 'miscalls', mode: 'copy'

        conda 'bedtools'

        input:
          val outputSuffix
          file clinvar_regions_for_methyl
		      file methyl_sorted
		      val n

        output:
          file "miscalls_${outputSuffix}_*.bed"

        shell:
        '''
        bedtools intersect -u -sorted -a !{methyl_sorted} -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' > miscalls_!{outputSuffix}_!{n.strip()}.bed
        '''
}

process intersectMiscalledGenes {

        conda 'bedtools'

        input:
                file methyl_miscalls
		file clinvar_annotations_low
		val n

        output:
                file 'miscalled_genes_dups_*.bed'
        shell:
        '''
        bedtools intersect -wa -sorted -a !{clinvar_annotations_low} -b !{methyl_miscalls} > miscalled_genes_dups_!{n.strip()}.bed
        '''
}

process countMiscalls {

        publishDir 'miscalls', mode: 'move'

        conda 'bedtools'

        input:
          val outputSuffix
                file miscalled_genes_dups_methyl
                val n

        output:
                file "miscalled_genes_${outputSuffix}_*.bed"
        shell:
        '''
        #!/usr/bin/env python
        import sys
        lines = {}
        out = open("miscalled_genes_!{outputSuffix}_!{n.strip()}.bed", "w")
        chr = ""


        def printLines(lines):
                for line in lines:
                        out.write(line+"\\t"+str(lines[line])+"\\n")

        with open("!{miscalled_genes_dups_methyl}") as f:
                for line in f:
                        if line.split()[0] != chr: #clear on new chrom to save memory and processing
                                printLines(lines)
                                lines = {}
                                chr = line.split()[0]
                        if line.strip() not in lines:
                                lines[line.strip()] = 0
                        lines[line.strip()] += 1
                        #print(line.strip())

                printLines(lines)
        '''
}

workflow {

  no_filter_mappability=params.no_filter_mappability

  bismap_bw = file(params.bismap)
  clinvar_regions = file(params.clinvarbed)
  ref_fai = file(params.ref_fai)
  ref = file(params.ref)


  basicPipeline(bismap_bbm, bismap_bw, ref, ref_fai, clinvar_regions, true, params.bams, params.suffix)


}

workflow basicPipeline {
  take:
    bismap_bbm
    bismap_bw
    ref
    ref_fai
    clinvar_regions
    no_filter_mappability
    pipeline_bams
    outputSuffix
  main:

  min_mapq = params.min_mapq
  bismap_cutoff = params.bismap_cutoff

  filter_opt = ""
  if(!no_filter_mappability)
  {
    filter_opt = "-B ${bismap_bbm}"
  }

  bam_names = pipeline_bams.toString().split(",")

  println(bam_names)

  bam_paths = bam_names.collect({ "readlink -f ${it}".execute().text })
  bam_list = Channel.from(bam_paths)

  bigWigToBedGraph(bismap_bw)
  convertFaiToGenome(ref_fai)
  complementBedgraph(convertFaiToGenome.out, bigWigToBedGraph.out)
  combineAllRegions(bigWigToBedGraph.out, complementBedgraph.out)
  getLowRegions(bigWigToBedGraph.out)
  combineLowRegions(getLowRegions.out, complementBedgraph.out)
  intersectMin(clinvar_regions, combineAllRegions.out)
  intersectMean(clinvar_regions, combineAllRegions.out)
  pasteMinMean(intersectMean.out, intersectMin.out)
  filterMinMean(pasteMinMean.out, bismap_cutoff)
  gsortMinMean(filterMinMean.out)
  intersectClinvarBismapLow(gsortMinMean.out, combineLowRegions.out)
  extractBam(ref, min_mapq, bam_list, filter_opt, ref_fai, bismap_bbm)
  stringParsing(extractBam.out[3])
  combineCalls(extractBam.out[0], extractBam.out[1], extractBam.out[2], stringParsing.out)
  intersectCalls(outputSuffix, intersectClinvarBismapLow.out, combineCalls.out, stringParsing.out)
  intersectMiscalledGenes(intersectCalls.out, gsortMinMean.out, stringParsing.out)
  countMiscalls(outputSuffix, intersectMiscalledGenes.out, stringParsing.out)
}
