nextflow.enable.dsl=2

params.clinvarBed="test_fixtures/chr18_clinvar_regions.bed"

params.bismap="test_fixtures/chr18_bismap.bw"

params.bismapBbm="test_fixtures/chr18_bismap.bbm"

params.refFai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.tmpdir = "/tmp"

params.minMapq=10

params.noFilterMappability=false

params.bismapCutoff=0.01

params.referenceName="ref"

params.aligner = "aligner"

params.seqMethod = "seq"

params.massMapping = "test_fixtures/mass_mapping.csv"

params.copy = false;

process bigWigToBedGraph {

    conda 'ucsc-bigwigtobedgraph=377'

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
        val tmpdir

    output:
        file 'ref.genome'

    shell:
    '''
    cut -f 1,2 !{ref_fai} | LC_ALL=C sort -k1,1 -T !{tmpdir}  > ref.genome
    '''
}

process extractBam {

    penv 'smp'
    conda 'methyldackel=0.6.0'
    cpus 8

    input:
        file ref
        val min_mapq
        each bam
        val filter_opt
        file ref_fai
        file bismap_bbm
  
    output:
        tuple file('*_CpG.bedGraph'), file('*_CHG.bedGraph'), file('*_CHH.bedGraph'), val("${bam.split("/")[-1].replace(".bam\n", "")}")
        
    shell:
    """
    MethylDackel extract -@ !{task.cpus} -q !{min_mapq} --CHH --CHG --OT 2,0,0,97 !{filter_opt} -o methyl_!{bam.split("/")[-1].replace(".bam\\n", "")} !{ref} !{bam.strip()}
    """
}

process combineCalls {

    conda 'coreutils=8.25'
    penv 'smp'
    cpus 8

    input:
        tuple file(methyl_cpg), file(methyl_chg), file(methyl_chh), val(n)
        val tmpdir

    output:
        tuple file('sorted_noheader_*.bedGraph'), val(n)
    
    shell:
    '''
    LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{methyl_cpg} !{methyl_chg} !{methyl_chh} | head -n -3 > sorted_noheader_!{n.strip()}.bedGraph
    '''
}

process complementBedgraph {

    conda 'bedtools=2.30.0'

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
    conda 'coreutils=8.25'

    input:
        file bismap_bedGraph_for_combine_all
        file zeroes_bedgraph
        val tmpdir

    output:
        file 'bismap_all.bedGraph'
    
    shell:
    '''
    LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{bismap_bedGraph_for_combine_all} !{zeroes_bedgraph} > bismap_all.bedGraph
    '''
}

process intersectMin {

    conda 'bedtools=2.30.0'

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

    conda 'bedtools=2.30.0'

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
    conda 'coreutils=8.25'

    input:
        file clinvar_regions_min_mean_map_low_unsorted
        val tmpdir

    output:
        file 'clinvar_annotations_low.bedGraph'
        
    shell:
    '''
    LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{clinvar_regions_min_mean_map_low_unsorted} > clinvar_annotations_low.bedGraph
    '''
}

process slopRegions {

    cpus 8
    penv 'smp'
    conda 'bedtools=2.30.0'

    input:
        file clinvar_regions_min_mean_map_low
        val tmpdir
        file ref_genome

    output:
        file 'clinvar_annotations_low_expanded.bedGraph'
        
    shell:
    '''
    cat !{clinvar_regions_min_mean_map_low} | awk '{ print $0"\tplaceholder" }' > tmp.bed; bedtools slop -l 2000 -r 0 -s -i tmp.bed -g !{ref_genome} | cut -f 1-12 > clinvar_annotations_low_expanded.bedGraph; rm tmp.bed
    '''
}


process getLowRegions {

    input:
        file bismap_bedGraph_for_awk_low_regions

    output:
        file 'bismap_low_nozeroes.bedGraph'
        
    shell:
    '''
    awk '$4 <= 0.01 { print }' !{bismap_bedGraph_for_awk_low_regions} > bismap_low_nozeroes.bedGraph
    '''
}



process combineLowRegions {

    cpus 8
    penv 'smp'
    conda 'coreutils=8.25'

    input:
        file low_bedgraph
        file zeroes_bedgraph
        val tmpdir

    output:
        file 'bismap_low.bedGraph'
        
    shell:
    '''
    LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{low_bedgraph} !{zeroes_bedgraph} > bismap_low.bedGraph
    '''
}

process intersectClinvarBismapLow {

    conda 'bedtools=2.30.0'

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

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        file regions_for_methyl
        tuple file(methyl_sorted), val(n)

    output:
      tuple file("lowmap_calls_${outputSuffix}_*.bed"), val(n)

    shell:
    '''
    bedtools intersect -u -sorted -a !{methyl_sorted} -b !{regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' > lowmap_calls_!{outputSuffix}_!{n.strip()}.bed
    '''
}

process annotateCalls {

    input:
        val outputSuffix
        tuple file(methyl_sorted), val(n)

    output:
      tuple file("annotated_calls_${outputSuffix}_*.bed"), val(n)

    shell:
    '''
    cat !{methyl_sorted} | awk '{ print $0"\t!{n.strip()}" }' > annotated_calls_!{outputSuffix}_!{n.strip()}.bed
    '''
}

process intersectClinvarCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        file clinvar_regions_for_methyl
        tuple file(methyl_sorted), val(n)

    output:
      tuple file("clinvar_calls_${outputSuffix}_*.bed"), val(n)

    shell:
    '''
    bedtools intersect -u -sorted -a !{methyl_sorted} -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' > clinvar_calls_!{outputSuffix}_!{n.strip()}.bed
    '''
}

process intersectClinvarCallsLowmap {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        file clinvar_regions_for_methyl
        tuple file(methyl_sorted), val(n)

    output:
      tuple file("clinvar_lowmap_calls_${outputSuffix}_*.bed"), val(n)

    shell:
    '''
    bedtools intersect -u -sorted -a !{methyl_sorted} -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' > clinvar_lowmap_calls_!{outputSuffix}_!{n.strip()}.bed
    '''
}


process intersectCalledGenes {

    conda 'bedtools=2.30.0'

    input:
        tuple file(methyl_calls), val(n)
        file clinvar_annotations_low

    output:
        tuple file('called_genes_dups_*.bed'), val(n)
        
    shell:
    '''
    bedtools intersect -wa -sorted -a !{clinvar_annotations_low} -b !{methyl_calls} > called_genes_dups_!{n.strip()}.bed
    '''
}


process intersectLowmapCalledGenes {

    conda 'bedtools=2.30.0'

    input:
        tuple file(methyl_lowmap_calls), val(n)
        file clinvar_annotations_low 

    output:
        tuple file('lowmap_called_genes_dups_*.bed'), val(n)
        
    shell:
    '''
    bedtools intersect -wa -sorted -a !{clinvar_annotations_low} -b !{methyl_lowmap_calls} > lowmap_called_genes_dups_!{n.strip()}.bed
    '''
}

process countLowmapGeneCalls {

    conda 'bedtools=2.30.0'
    //use nextflow bin dir for script?
    input:
        val outputSuffix
        tuple file(lowmap_called_genes_dups_methyl), val(n)
        val basePath

    output:
        tuple file("lowmap_called_genes_${outputSuffix}_*.bed"), val(n)
        
    shell:
    '''
    python !{basePath}/count_calls.py lowmap_called_genes_!{outputSuffix}_!{n.strip()}.bed !{lowmap_called_genes_dups_methyl} !{n.strip()}
    '''
}

process countGeneCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple file(called_genes_dups_methyl), val(n)
        val basePath

    output:
        tuple file("called_genes_${outputSuffix}_*.bed"), val(n)
        
    shell:
    '''
    python !{basePath}/count_calls.py called_genes_!{outputSuffix}_!{n.strip()}.bed !{called_genes_dups_methyl} !{n.strip()}
    '''
}


process addParamsToLowmapCalledGenes {

    publishDir 'methyl_calls/lowmap_called_genes', mode: "${publish_mode}", pattern: "lowmap_called_genes_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(lowmap_called_genes), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("lowmap_called_genes_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{lowmap_called_genes} lowmap_called_genes_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 13
    '''
}

process addParamsToCalledGenes {

    publishDir 'methyl_calls/called_genes', mode: "${publish_mode}", pattern: "called_genes_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(called_genes), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("called_genes_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{called_genes} called_genes_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 13
    '''
}

process addParamsToCalls {

    publishDir 'methyl_calls/calls', mode: "${publish_mode}", pattern: "calls_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(calls), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("calls_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{calls} calls_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process addParamsToLowmapCalls {

    publishDir 'methyl_calls/lowmap_calls', mode: "${publish_mode}", pattern: "lowmap_calls_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(lowmap_calls), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("lowmap_calls_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{lowmap_calls} lowmap_calls_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process addParamsToClinvarCalls {

    publishDir 'methyl_calls/clinvar_calls', mode: "${publish_mode}", pattern: "calls_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(lowmap_calls), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("calls_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{lowmap_calls} calls_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process addParamsToClinvarLowmapCalls {

    publishDir 'methyl_calls/clinvar_lowmap_calls', mode: "${publish_mode}", pattern: "lowmap_calls_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(lowmap_calls), val(n)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(n), file("lowmap_calls_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{lowmap_calls} lowmap_calls_!{outputSuffix}_!{n.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process getPureMethylSites {

    publishDir 'methyl_calls/pure_mixed_calls', mode: "${publish_mode}", pattern: "pure_methyl_calls_${outputSuffix}_*.tsv"
    conda 'bedtools=2.30.0'

    input:
       tuple val(n), file(calls)
       val outputSuffix
       val publish_mode
    output:
        tuple val(n), file("pure_methyl_calls_${outputSuffix}_*.tsv"), file(calls)
    shell:
    '''
    awk '($4 == 100 || $4 == 0) { print }' !{calls} > pure_methyl_calls_!{outputSuffix}_!{n.strip()}.tsv
    '''
}

process getMixedMethylSites {

    cpus 2
    penv 'smp'
    publishDir 'methyl_calls/pure_mixed_calls', mode: "${publish_mode}", pattern: "mixed_methyl_calls_${outputSuffix}_*.tsv"
    conda 'bedtools=2.30.0'

    input:
       tuple val(n), file(pure_calls), file(calls)
       val outputSuffix
       val publish_mode
    output:
        tuple val(n), file("mixed_methyl_calls_${outputSuffix}_*.tsv")
    shell:
    '''
    bedtools subtract -A -a !{calls} -b !{pure_calls} > mixed_methyl_calls_!{outputSuffix}_!{n.strip()}.tsv
    '''
}

workflow {

    no_filter_mappability=params.noFilterMappability

    bismap_bw = file(params.bismap)
    clinvar_regions = file(params.clinvarBed)
    ref_fai = file(params.refFai)
    ref = file(params.ref)
    bismap_bbm= file(params.bismapBbm)
    mass_mapping = file(params.massMapping)

    publishMode = "symlink"

    if(params.copy)
    {
        publishMode = "copy"
    }

    refPreproc(bismap_bw, ref, ref_fai, clinvar_regions, params.tmpdir, params.minMapq, params.bismapCutoff)
    
    basicPipeline(bismap_bbm, ref, ref_fai, no_filter_mappability, params.bams, params.referenceName, params.aligner, params.seqMethod, params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, refPreproc.out[0], refPreproc.out[1], refPreproc.out[2], publishMode)


}

workflow refPreproc {
    take:
        bismap_bw
        ref
        ref_fai
        clinvar_regions
        tmpdir
        min_mapq
        bismap_cutoff
    main:
        bigWigToBedGraph(bismap_bw)
        convertFaiToGenome(ref_fai, tmpdir)
        complementBedgraph(convertFaiToGenome.out, bigWigToBedGraph.out)
        combineAllRegions(bigWigToBedGraph.out, complementBedgraph.out, tmpdir)
        getLowRegions(bigWigToBedGraph.out)
        combineLowRegions(getLowRegions.out, complementBedgraph.out, tmpdir)
        intersectMin(clinvar_regions, combineAllRegions.out)
        intersectMean(clinvar_regions, combineAllRegions.out)
        pasteMinMean(intersectMean.out, intersectMin.out)
        filterMinMean(pasteMinMean.out, bismap_cutoff)
        gsortMinMean(filterMinMean.out, tmpdir)
        slopRegions(gsortMinMean.out, tmpdir, convertFaiToGenome.out)
        intersectClinvarBismapLow(slopRegions.out, combineLowRegions.out)
    emit:
        intersectClinvarBismapLow.out
        combineLowRegions.out
        slopRegions.out
}

workflow basicPipeline {
    take:
        bismap_bbm
        ref
        ref_fai
        no_filter_mappability
        pipeline_bams
        referenceName
        aligner
        seqMethod
        tmpdir
        min_mapq
        bismap_cutoff
        massMapping
        clinvarBismapLow
        lowRegions
        lowGenes
        publish_mode
    main:

        filterStr = "nofilt"

        filter_opt = ""
        
        if(!no_filter_mappability)
        {
            filter_opt = "-B ${bismap_bbm}"
            filterStr = "filt"
        }

        outputSuffix = referenceName+"_"+aligner+"_"+seqMethod+"_"+filterStr

        bam_names = pipeline_bams.toString().split(",")

        println(bam_names)

        bam_paths = bam_names.collect({ "readlink -f ${it}".execute().text })
        bam_list = Channel.from(bam_paths)

        extractBam(ref, min_mapq, bam_list, filter_opt, ref_fai, bismap_bbm)
        combineCalls(extractBam.out, tmpdir)
        intersectClinvarCalls(outputSuffix, lowGenes, combineCalls.out)
        intersectClinvarCallsLowmap(outputSuffix, clinvarBismapLow, combineCalls.out)
        intersectCalls(outputSuffix, lowRegions, combineCalls.out)
        annotateCalls(outputSuffix, combineCalls.out)
        intersectLowmapCalledGenes(intersectCalls.out, lowGenes)
        intersectCalledGenes(annotateCalls.out, lowGenes)
        countLowmapGeneCalls(outputSuffix, intersectLowmapCalledGenes.out, workflow.projectDir)
        countGeneCalls(outputSuffix, intersectCalledGenes.out, workflow.projectDir)
        addParamsToCalls(outputSuffix, annotateCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToLowmapCalls(outputSuffix, intersectCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToClinvarCalls(outputSuffix, intersectClinvarCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToClinvarLowmapCalls(outputSuffix, intersectClinvarCallsLowmap.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToCalledGenes(outputSuffix, countGeneCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToLowmapCalledGenes(outputSuffix, countLowmapGeneCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        getPureMethylSites(addParamsToCalls.out, outputSuffix, publish_mode)
        getMixedMethylSites(getPureMethylSites.out, outputSuffix, publish_mode)
    emit:
        addParamsToLowmapCalls.out
        addParamsToCalls.out
        addParamsToClinvarLowmapCalls.out
        getPureMethylSites.out
        getMixedMethylSites.out
        addParamsToClinvarCalls.out
}
