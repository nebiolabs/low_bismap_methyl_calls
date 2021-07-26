nextflow.enable.dsl=2

include { refPreproc as PreprocRef } from './miscalls_one.nf'

include { basicPipeline as NoFilt } from './miscalls_one.nf'
include { basicPipeline as Filt } from './miscalls_one.nf'

params.clinvarBed="test_fixtures/chr18_clinvar_regions.bed"

params.bismapBw="test_fixtures/chr18_bismap.bw"

params.bismapBbm="test_fixtures/chr18_bismap.bbm"

params.refFai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.minMapq=10

params.bismapCutoff=0.01

params.tmpdir = "/tmp"

params.massMapping = "test_fixtures/mass_mapping.csv"

params.referenceName="ref"

params.aligner = "aligner"

params.seqMethod = "seq"

params.copy = false

process combineFilteringLowmapCalls {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "miscalls_lowmap_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(nofilt), file(filt)
        val publish_mode

    output:
        tuple val(groupKey), file("miscalls_lowmap_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -v -a !{nofilt} -b !{filt} | cut -f 1-10,12 > miscalls_lowmap_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineFilteringClinvarLowmapCalls {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "clinvar_miscalls_lowmap_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(nofilt), file(filt)
        val publish_mode

    output:
        tuple val(groupKey), file("clinvar_miscalls_lowmap_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -v -a !{nofilt} -b !{filt} | cut -f 1-10,12 > clinvar_miscalls_lowmap_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineFilteringClinvarCalls {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "clinvar_miscalls_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(nofilt), file(filt)
        val publish_mode

    output:
        tuple val(groupKey), file("clinvar_miscalls_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -v -a !{nofilt} -b !{filt} | cut -f 1-10,12 > clinvar_miscalls_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineFilteringCalls {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "miscalls_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(nofilt), file(filt)
        val publish_mode

    output:
        tuple val(groupKey), file("miscalls_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -v -a !{nofilt} -b !{filt} | cut -f 1-10,12 > miscalls_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineFilteringClinvarGenesCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(calls)
        file clinvar_data

    output:
        tuple val(groupKey), file("miscalled_clinvar_genes_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -wa -sorted -a !{clinvar_data} -b !{calls} > miscalled_clinvar_genes_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineFilteringClinvarGenesLowmapCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(lowmap_miscalls)
        file clinvar_data

    output:
        tuple val(groupKey), file("miscalled_clinvar_genes_lowmap_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -wa -sorted -a !{clinvar_data} -b !{lowmap_miscalls} > miscalled_clinvar_genes_lowmap_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process combineClinvarGenesResolvedCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(resolved_calls)
        file clinvar_data

    output:
        tuple val(groupKey), file("resolved_clinvar_genes_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -wa -sorted -a !{clinvar_data} -b !{resolved_calls} > resolved_clinvar_genes_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process countClinvarGeneResolvedCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(resolved_genes_dups_methyl)
        val basePath

    output:
        tuple file("resolved_clinvar_genes_${outputSuffix}_*.bed"), val(groupKey)
        
    shell:
    '''
    python !{basePath}/count_calls.py resolved_clinvar_genes_!{outputSuffix}_!{groupKey}.bed !{resolved_genes_dups_methyl} !{groupKey.strip()}
    '''
}


process countClinvarGeneLowmapCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(lowmap_miscalled_genes_dups_methyl)
        val basePath

    output:
        tuple file("miscalled_genes_lowmap_${outputSuffix}_*.bed"), val(groupKey)
        
    shell:
    '''
    python !{basePath}/count_calls.py miscalled_genes_lowmap_!{outputSuffix}_!{groupKey}.bed !{lowmap_miscalled_genes_dups_methyl} !{groupKey.strip()}
    '''
}

process countClinvarGeneCalls {

    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(miscalled_genes_dups_methyl)
        val basePath

    output:
        tuple file("miscalled_genes_${outputSuffix}_*.bed"), val(groupKey)
        
    shell:
    '''
    python !{basePath}/count_calls.py miscalled_genes_!{outputSuffix}_!{groupKey}.bed !{miscalled_genes_dups_methyl} !{groupKey.strip()}
    '''
}

process addParamsToClinvarLowmapCallsCombined {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "miscalled_clinvar_genes_lowmap_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(lowmap_miscalls), val(groupKey)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(groupKey), file("miscalled_clinvar_genes_lowmap_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{lowmap_miscalls} miscalled_clinvar_genes_lowmap_!{outputSuffix}_!{groupKey.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process addParamsToClinvarCallsCombined {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "miscalled_clinvar_genes_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(miscalls), val(groupKey)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(groupKey), file("miscalled_clinvar_genes_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{miscalls} miscalled_clinvar_genes_!{outputSuffix}_!{groupKey.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process addParamsToClinvarResolvedCallsCombined {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "resolved_clinvar_genes_${outputSuffix}_*.tsv"

    input:
        val outputSuffix
        tuple file(resolved_calls), val(groupKey)
        val referenceName
        val aligner
        val seqMethod
        val filterStr
        file massMapping
        val basePath
        val publish_mode

    output:
        tuple val(groupKey), file("resolved_clinvar_genes_${outputSuffix}_*.tsv")
        
    shell:
    '''
    python !{basePath}/add_params.py !{resolved_calls} resolved_clinvar_genes_!{outputSuffix}_!{groupKey.strip()}.tsv !{massMapping} !{referenceName} !{aligner} !{seqMethod} !{filterStr} 6
    '''
}

process combineResolvedCalls {

    publishDir 'methyl_calls/miscalls', mode: "${publish_mode}", pattern: "resolved_calls_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(nofilt), file(filt)
        val publish_mode

    output:
        tuple val(groupKey), file("resolved_calls_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -wa -wb -a !{nofilt} -b !{filt} | cut -f 1-6,16-24 > resolved_calls_!{outputSuffix}_!{groupKey}.tsv
    '''
}

process calcUnderfilteredCalls {

    publishDir 'methyl_calls/underfiltered_calls', mode: "${publish_mode}", pattern: "underfiltered_calls_${outputSuffix}_${groupKey}.tsv"
    conda 'bedtools=2.30.0'

    input:
        val outputSuffix
        tuple val(groupKey), file(resolved_calls)
        file lowmap_data
        val publish_mode

    output:
        tuple val(groupKey), file("underfiltered_calls_${outputSuffix}_${groupKey}.tsv")

    shell:
    '''
    bedtools intersect -sorted -v -a !{resolved_calls} -b !{lowmap_data} > underfiltered_calls_!{outputSuffix}_!{groupKey}.tsv
    '''
}

workflow {

    bismap_bw = file(params.bismapBw)
    bismap_bbm= file(params.bismapBbm)
    clinvar_regions = file(params.clinvarBed)
    ref_fai = file(params.refFai)
    ref = file(params.ref)

    mass_mapping = file(params.massMapping)

    publishMode = "symlink"

    if(params.copy)
    {
        publishMode = "copy"
    }

    PreprocRef(bismap_bw, ref, ref_fai, clinvar_regions, params.tmpdir, params.minMapq, params.bismapCutoff)
    runPair(bismap_bbm, ref, ref_fai, params.bams, params.referenceName, params.aligner, params.seqMethod, params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRef.out[0], PreprocRef.out[1], PreprocRef.out[2], publishMode)

}

workflow runPair {

    take:
        bismap_bbm
        ref
        ref_fai
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

        outputSuffix = referenceName+"_"+aligner+"_"+seqMethod

        filterStr = "combined"

        NoFilt(bismap_bbm, ref, ref_fai, true, pipeline_bams, referenceName, aligner, seqMethod, tmpdir, min_mapq, bismap_cutoff, massMapping, clinvarBismapLow, lowRegions, lowGenes, publish_mode)
        Filt(bismap_bbm, ref, ref_fai, false, pipeline_bams, referenceName, aligner, seqMethod, tmpdir, min_mapq, bismap_cutoff, massMapping, clinvarBismapLow, lowRegions, lowGenes, publish_mode)
        lowmapCallsChannel = NoFilt.out[0].join(Filt.out[0])
        callsChannel = NoFilt.out[1].join(Filt.out[1])
        clinvarLowmapCallsChannel = NoFilt.out[2].join(Filt.out[2])
        clinvarCallsChannel = NoFilt.out[5].join(Filt.out[5])
        callsMixedPureChannel = NoFilt.out[4].join(Filt.out[3])

        combineFilteringLowmapCalls(outputSuffix, lowmapCallsChannel, publish_mode)
        combineFilteringCalls(outputSuffix, callsChannel, publish_mode)
        combineFilteringClinvarLowmapCalls(outputSuffix, clinvarLowmapCallsChannel, publish_mode)
        combineFilteringClinvarCalls(outputSuffix, clinvarCallsChannel, publish_mode)

        calcUnderfilteredCalls(outputSuffix, combineFilteringCalls.out, clinvarBismapLow, publish_mode)

        combineFilteringClinvarGenesCalls(outputSuffix, combineFilteringCalls.out, lowGenes)
        combineFilteringClinvarGenesLowmapCalls(outputSuffix, combineFilteringLowmapCalls.out, lowGenes)


        countClinvarGeneCalls(outputSuffix, combineFilteringClinvarGenesCalls.out, workflow.projectDir)
        countClinvarGeneLowmapCalls(outputSuffix, combineFilteringClinvarGenesLowmapCalls.out, workflow.projectDir)


        addParamsToClinvarCallsCombined(outputSuffix, countClinvarGeneCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
        addParamsToClinvarLowmapCallsCombined(outputSuffix, countClinvarGeneLowmapCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)

        combineResolvedCalls(outputSuffix, callsMixedPureChannel, publish_mode)

        combineClinvarGenesResolvedCalls(outputSuffix, combineResolvedCalls.out, lowGenes)
        countClinvarGeneResolvedCalls(outputSuffix, combineClinvarGenesResolvedCalls.out, workflow.projectDir)
        addParamsToClinvarResolvedCallsCombined(outputSuffix, countClinvarGeneResolvedCalls.out, referenceName, aligner, seqMethod, filterStr, massMapping, workflow.projectDir, publish_mode)
}