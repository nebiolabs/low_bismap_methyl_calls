params.clinvar="test_fixtures/chr18_clinvar.vcf"

params.gencode="test_fixtures/chr18_gencode.gff"

params.bismap="test_fixtures/chr18_bismap.bw"

params.ref_fai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams="test_fixtures/*.bam"

params.min_mapq=20

params.no_filter_mappability=false

no_filter_mappability=params.no_filter_mappability

params.bismap_cutoff=0.01

sortmem = params.sortmem

clinvar_vcf = file(params.clinvar)

gencode_gff = file(params.gencode)

bismap_bw = file(params.bismap)

bismap_bw_for_methyldackel = file(params.bismap)

filter_opt = ""

if(!no_filter_mappability)
{
	filter_opt = "-M ${bismap_bw_for_methyldackel}"
}

ref_fai = file(params.ref_fai)

ref = file(params.ref)

min_mapq = params.min_mapq

bismap_cutoff = params.bismap_cutoff

bam_names = params.bams.toString().split(",")

println(bam_names)

f_name = Channel.create()

println("[DEBUG] 1")

bam_paths = bam_names.collect({ "readlink -f ${it}".execute().text })

println("[DEBUG] 2")

bam_list = Channel.from(bam_paths)

println("[DEBUG] About to start!")
process vcf2bed {
	
	input:
		file clinvar_vcf

	output:
		file 'clinvar_bed.bed' into clinvar_bed

	shell:
	'''
	vcf2bed < !{clinvar_vcf} | sed "s/^MT/M/g" | sed -E "s/(^[0-9]+)/chr\\1/g" > clinvar_bed.bed
	'''
}


process gff2bed {

	input:
		file gencode_gff

	output:
		file 'gencode_bed.bed' into gencode_bed

	shell:
	'''
	gff2bed < !{gencode_gff} > gencode_bed.bed
	'''
}

process bigWigToBedGraph {

        input:
                file bismap_bw

        output:
                file 'bismap.bedGraph' into bismap_bedGraph_for_combine_all
		file 'bismap.bedGraph' into bismap_bedGraph_for_complement
		file 'bismap.bedGraph' into bismap_bedGraph_for_awk_low_regions

        shell:
        '''
        bigWigToBedGraph !{bismap_bw} bismap.bedGraph
        '''
}

process convertFaiToGenome {

        input:
                file ref_fai

        output:
                file 'ref.genome' into ref_genome

        shell:
        '''
        cut -f 1,2 !{ref_fai} > ref.genome
        '''
}

process extractBam {

	cpus 8

        input:
		val min_mapq
                each bam from bam_list
		val filter_opt
	
        output:
                file '*_CpG.bedGraph' into methyl_cpg
		file '*_CHG.bedGraph' into methyl_chg
		file '*_CHH.bedGraph' into methyl_chh
		val bam into f_name_raw
        shell:
        """
        MethylDackel extract -@ !{task.cpus} -q !{min_mapq} --CHH --CHG --OT 2,0,0,97 !{filter_opt} -o methyl_!{bam.split("/")[-1].replace(".bam\\n", "")} !{ref} !{bam} 
        """
}

process stringParsing {

        input:
		val bam_raw from f_name_raw

	exec:
	f_name << "$bam_raw".split("/")[-1].replace(".bam\n", "")
}



process combineCalls {

	cpus 8

        input:
                file methyl_cpg
		file methyl_chg
		file methyl_chh
		val n from f_name

        output:
                file 'sorted_noheader_*.bedGraph' into methyl_sorted
		val n into f_name_2
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{methyl_cpg} !{methyl_chg} !{methyl_chh} | head -n -3 > sorted_noheader_!{n}.bedGraph
        '''
}

process intersectClinvarGencode {

        input:
                file clinvar_bed
		file gencode_bed

        output:
                file 'clinvar_regions.bed' into clinvar_regions_for_min
		file 'clinvar_regions.bed' into clinvar_regions_for_mean
        shell:
        '''
        bedtools intersect -u -a !{gencode_bed} -b !{clinvar_bed} > clinvar_regions.bed
        '''
}

process complementBedgraph {

        input:
                file ref_genome
		file bismap_bedGraph_for_complement

        output:
                file 'bismap_zeroes.bedGraph' into zeroes_bedgraph
        shell:
        '''
        bedtools complement -i !{bismap_bedGraph_for_complement} -g !{ref_genome} | awk '{ print $0"\t"0 }'> bismap_zeroes.bedGraph
        '''
}

process combineAllRegions {

	cpus 8

        input:

		file bismap_bedGraph_for_combine_all
		file zeroes_bedgraph

        output:
                file 'bismap_all.bedGraph' into bismap_bedgraph_for_mean, bismap_bedgraph_for_min
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{bismap_bedGraph_for_combine_all} !{zeroes_bedgraph} > bismap_all.bedGraph
        '''
}

process intersectMin {

        input:
                file clinvar_regions_for_min
		file bismap_bedgraph_for_min

        output:
                file 'clinvar_regions_min.bedGraph' into clinvar_regions_min_map
        shell:
        '''
        bedtools map -o min -c 4 -a !{clinvar_regions_for_min} -b !{bismap_bedgraph_for_min} > clinvar_regions_min.bedGraph
        '''
}

process intersectMean {

        input:
                file clinvar_regions_for_mean
		file bismap_bedgraph_for_mean

        output:
                file 'clinvar_regions_mean.bedGraph' into clinvar_regions_mean_map
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
                file 'clinvar_regions_min_mean.bedGraph' into clinvar_regions_min_mean_map
        shell:
        '''
        paste !{clinvar_regions_min_map} <(cut -f 11 !{clinvar_regions_mean_map}) > clinvar_regions_min_mean.bedGraph
        '''
}

process filterMinMean {

        input:
                file clinvar_regions_min_mean_map

        output:
                file 'clinvar_regions_min_mean_low.bedGraph' into clinvar_regions_min_mean_map_low_unsorted
        shell:
        '''
        awk '$11 <= !{bismap_cutoff} { print }' !{clinvar_regions_min_mean_map} > clinvar_regions_min_mean_low.bedGraph
        '''
}

process gsortMinMean {

	cpus 8

        input:
                file clinvar_regions_min_mean_map_low_unsorted

        output:
                file 'clinvar_annotations_low.bedGraph' into clinvar_annotations_low
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{clinvar_regions_min_mean_map_low_unsorted} > clinvar_annotations_low.bedGraph
        '''
}

process getLowRegions {

        input:

		file bismap_bedGraph_for_awk_low_regions

        output:
                file 'bismap_low.bedGraph' into low_bedgraph
        shell:
        '''
        awk '$4 <= 0.01 { print }' !{bismap_bedGraph_for_awk_low_regions} > bismap_low.bedGraph
        '''
}



process combineLowRegions {

	cpus 8

        input:

		file low_bedgraph
		file zeroes_bedgraph

        output:
                file 'bismap_low.bedGraph' into bismap_low
        shell:
        '''
        sort -k1,1 -k2,2n -S 32G --parallel !{task.cpus} !{low_bedgraph} !{zeroes_bedgraph} > bismap_low.bedGraph
        '''
}

process intersectClinvarBismapLow {

        input:
                file clinvar_annotations_low
		file bismap_low

        output:
		file 'clinvar_low_regions.bed' into clinvar_regions_for_methyl
        shell:
        '''
        bedtools intersect -sorted -a !{bismap_low} -b !{clinvar_annotations_low} > clinvar_low_regions.bed
        '''
}

process intersectCalls {

	publishDir 'miscalls', mode: 'copy'

        input:
                file clinvar_regions_for_methyl
		file methyl_sorted
		val n from f_name_2

        output:
                file 'miscalls_*.bed' into methyl_miscalls
		val n into f_name_3
        shell:
        '''
        bedtools intersect -u -sorted -a !{methyl_sorted} -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n}" }' > miscalls_!{n}.bed
        '''
}

process intersectMiscalledGenes {

        input:
                file methyl_miscalls
		file clinvar_annotations_low
		val n from f_name_3

        output:
                file 'miscalled_genes_dups_*.bed' into miscalled_genes_dups_methyl
		val n into f_name_4
        shell:
        '''
        bedtools intersect -wa -sorted -a !{clinvar_annotations_low} -b !{methyl_miscalls} > miscalled_genes_dups_!{n}.bed
        '''
}

process countMiscalls {

        publishDir 'miscalls', mode: 'move'

        input:
                file miscalled_genes_dups_methyl
                val n from f_name_4

        output:
                file 'miscalled_genes_*.bed' into miscalled_genes_methyl
        shell:
        '''
        #!/usr/bin/env python
        import sys
        lines = {}
        out = open("miscalled_genes_!{n}.bed", "w")
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

