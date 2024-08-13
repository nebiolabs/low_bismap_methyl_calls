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

params.hardlink = false

params.outputDir="methyl_calls"

outputDir=params.outputDir

process combineFilteringLowmapCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalls", mode: "${publish_mode}", pattern: "miscalls_lowmap_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 coreutils=9.1 zstd=1.5.6'

        input:
                tuple val(groupKey), file(nofilt), file(filt)
                val publish_mode

        output:
                tuple val(groupKey), file("miscalls_lowmap_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -v -a <(zstd -d -f --stdout !{nofilt}) -b <(zstd -d -f --stdout !{filt}) | cut -f 1-10,12 | zstd -o miscalls_lowmap_!{groupKey}.tsv.zst
                '''
}

process combineFilteringClinvarLowmapCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalls", mode: "${publish_mode}", pattern: "clinvar_miscalls_lowmap_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 coreutils=9.1 zstd=1.5.6'

        input:
                tuple val(groupKey), file(nofilt), file(filt)
                val publish_mode

        output:
                tuple val(groupKey), file("clinvar_miscalls_lowmap_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -v -a <(zstd -d -f --stdout !{nofilt}) -b <(zstd -d -f --stdout !{filt}) | cut -f 1-10,12 | zstd -o clinvar_miscalls_lowmap_!{groupKey}.tsv.zst
                '''
}

process combineFilteringClinvarCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalls", mode: "${publish_mode}", pattern: "clinvar_miscalls_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 coreutils=9.1 zstd=1.5.6'

        input:
                tuple val(groupKey), file(nofilt), file(filt)
                val publish_mode

        output:
                tuple val(groupKey), file("clinvar_miscalls_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -v -a <(zstd -d -f --stdout !{nofilt}) -b <(zstd -d -f --stdout !{filt}) | cut -f 1-10,12 | zstd -o clinvar_miscalls_!{groupKey}.tsv.zst
                '''
}

process combineFilteringCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalls", mode: "${publish_mode}", pattern: "miscalls_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 coreutils=9.1 zstd=1.5.6'

        input:
                tuple val(groupKey), file(nofilt), file(filt)
                val publish_mode

        output:
                tuple val(groupKey), file("miscalls_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -v -a <(zstd -d -f --stdout !{nofilt}) -b <(zstd -d -f --stdout !{filt}) | cut -f 1-10,12 | zstd -o miscalls_!{groupKey}.tsv.zst
                '''
}

process combineFilteringClinvarGenesCalls {
        
        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple val(groupKey), file(calls)
                file clinvar_data

        output:
                tuple val(groupKey), file("miscalled_clinvar_genes_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -wa -sorted -a !{clinvar_data} -b <(zstd -d -f --stdout !{calls}) | zstd -o miscalled_clinvar_genes_!{groupKey}.tsv.zst
                '''
}

process combineFilteringClinvarGenesLowmapCalls {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple val(groupKey), file(lowmap_miscalls)
                file clinvar_data

        output:
                tuple val(groupKey), file("miscalled_clinvar_genes_lowmap_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -wa -sorted -a !{clinvar_data} -b <(zstd -d -f --stdout !{lowmap_miscalls}) | zstd -o miscalled_clinvar_genes_lowmap_!{groupKey}.tsv.zst
                '''
}

process combineClinvarGenesResolvedCalls {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple val(groupKey), file(resolved_calls)
                file clinvar_data

        output:
                tuple val(groupKey), file("resolved_clinvar_genes_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -wa -sorted -a !{clinvar_data} -b <(zstd -d -f --stdout !{resolved_calls}) | zstd -o resolved_clinvar_genes_!{groupKey}.tsv.zst
                '''
}

process countClinvarGeneResolvedCalls {

        penv 'smp'
        cpus 1
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple val(groupKey), file(resolved_genes_dups_methyl)
                val basePath

        output:
                tuple file("resolved_clinvar_genes_*.bed.zst"), val(groupKey)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{resolved_genes_dups_methyl} | python !{basePath}/count_calls.py !{groupKey.strip()} | zstd -o resolved_clinvar_genes_!{groupKey}.bed.zst
                '''
}


process countClinvarGeneLowmapCalls {

        penv 'smp'
        cpus 1
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple val(groupKey), file(lowmap_miscalled_genes_dups_methyl)
                val basePath

        output:
                tuple file("miscalled_genes_lowmap_*.bed.zst"), val(groupKey)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_miscalled_genes_dups_methyl} | python !{basePath}/count_calls.py !{groupKey.strip()} | zstd -o miscalled_genes_lowmap_!{groupKey}.bed.zst
                '''
}

process countClinvarGeneCalls {

        penv 'smp'
        cpus 1
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple val(groupKey), file(miscalled_genes_dups_methyl)
                val basePath

        output:
                tuple file("miscalled_genes_*.bed.zst"), val(groupKey)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{miscalled_genes_dups_methyl} | python !{basePath}/count_calls.py !{groupKey.strip()} | zstd -o miscalled_genes_!{groupKey}.bed.zst 
                '''
}

process addParamsToClinvarLowmapCallsCombined {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalled_genes", mode: "${publish_mode}", pattern: "miscalled_clinvar_genes_lowmap_*.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(lowmap_miscalls), val(groupKey)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(groupKey), file("miscalled_clinvar_genes_lowmap_*.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_miscalls} | python !{basePath}/add_params_new.py !{filterStr} 13 | zstd -o miscalled_clinvar_genes_lowmap_!{groupKey.strip()}.tsv.zst
                '''
}

process addParamsToClinvarCallsCombined {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalled_genes", mode: "${publish_mode}", pattern: "miscalled_clinvar_genes_*.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(miscalls), val(groupKey)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(groupKey), file("miscalled_clinvar_genes_*.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{miscalls} | python !{basePath}/add_params_new.py !{filterStr} 13 | zstd -o miscalled_clinvar_genes_!{groupKey.strip()}.tsv.zst 
                '''
}

process addParamsToClinvarResolvedCallsCombined {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalled_genes", mode: "${publish_mode}", pattern: "resolved_clinvar_genes_*.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(resolved_calls), val(groupKey)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(groupKey), file("resolved_clinvar_genes_*.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{resolved_calls} | python !{basePath}/add_params_new.py !{filterStr} 13 | zstd -o resolved_clinvar_genes_!{groupKey.strip()}.tsv.zst  
                '''
}

process combineResolvedCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/miscalls", mode: "${publish_mode}", pattern: "resolved_calls_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 coreutils=9.1 zstd=1.5.6'

        input:
                tuple val(groupKey), file(nofilt), file(filt)
                val publish_mode

        output:
                tuple val(groupKey), file("resolved_calls_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -wa -wb -a <(zstd -d -f --stdout !{nofilt}) -b <(zstd -d -f --stdout !{filt}) | cut -f 1-6,16-24 | zstd -o resolved_calls_!{groupKey}.tsv.zst
                '''
}

process calcUnderfilteredCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/underfiltered_calls", mode: "${publish_mode}", pattern: "underfiltered_calls_${groupKey}.tsv.zst"
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple val(groupKey), file(resolved_calls)
                file lowmap_data
                val publish_mode

        output:
                tuple val(groupKey), file("underfiltered_calls_${groupKey}.tsv.zst")

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -v -a <(zstd -d -f --stdout !{resolved_calls}) -b !{lowmap_data} | zstd -o underfiltered_calls_!{groupKey}.tsv.zst
                '''
}

workflow {

    bismap_bw = file(params.bismapBw)
    bismap_bbm= file(params.bismapBbm)
    clinvar_regions = file(params.clinvarBed)
    ref_fai = file(params.refFai)
    ref = file(params.ref)

    pipeline_bams = Channel.fromFilePairs(params.bam_glob+"{.bam,.bam.bai}")

    publishMode = "symlink"

    if(params.hardlink)
    {
        publishMode = "link"
    }

    PreprocRef(bismap_bw, ref, ref_fai, clinvar_regions, params.tmpdir, params.minMapq, params.bismapCutoff)
    runPair(bismap_bbm, ref, ref_fai, pipeline_bams, params.tmpdir, params.minMapq, params.bismapCutoff, PreprocRef.out[0], PreprocRef.out[1], PreprocRef.out[2], publishMode)

}

workflow runPair {

    take:
        bismap_bbm
        ref
        ref_fai
        pipeline_bams
        tmpdir
        min_mapq
        bismap_cutoff
        clinvarBismapLow
        lowRegions
        lowGenes
        publish_mode
    main:

        //seqParams = pipeline_seqname.split("\\.")[1..-1]

        //seqName = seqParams[0]
        //inputMass = seqParams[1]
        //seqMethod = seqParams[2]
        //referenceName = seqParams[3]
        //aligner = seqParams[4]

        //outputSuffix = seqName+"_"+referenceName+"_"+aligner

        filterStr = "combined"

        NoFilt(bismap_bbm, ref, ref_fai, true, pipeline_bams, tmpdir, min_mapq, bismap_cutoff, clinvarBismapLow, lowRegions, lowGenes, publish_mode)
        Filt(bismap_bbm, ref, ref_fai, false, pipeline_bams, tmpdir, min_mapq, bismap_cutoff, clinvarBismapLow, lowRegions, lowGenes, publish_mode)
        lowmapCallsChannel = NoFilt.out[0].join(Filt.out[0])
        callsChannel = NoFilt.out[1].join(Filt.out[1])
        clinvarLowmapCallsChannel = NoFilt.out[2].join(Filt.out[2])
        clinvarCallsChannel = NoFilt.out[5].join(Filt.out[5])
        callsMixedPureChannel = NoFilt.out[4].join(Filt.out[3])

        combineFilteringLowmapCalls(lowmapCallsChannel, publish_mode)
        combineFilteringCalls(callsChannel, publish_mode)
        combineFilteringClinvarLowmapCalls(clinvarLowmapCallsChannel, publish_mode)
        combineFilteringClinvarCalls(clinvarCallsChannel, publish_mode)

        calcUnderfilteredCalls(combineFilteringCalls.out, lowRegions, publish_mode)

        combineFilteringClinvarGenesCalls(combineFilteringCalls.out, lowGenes)
        combineFilteringClinvarGenesLowmapCalls(combineFilteringLowmapCalls.out, lowGenes)


        countClinvarGeneCalls(combineFilteringClinvarGenesCalls.out, workflow.projectDir)
        countClinvarGeneLowmapCalls(combineFilteringClinvarGenesLowmapCalls.out, workflow.projectDir)


        addParamsToClinvarCallsCombined(countClinvarGeneCalls.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToClinvarLowmapCallsCombined(countClinvarGeneLowmapCalls.out, filterStr, workflow.projectDir, publish_mode)

        combineResolvedCalls(callsMixedPureChannel, publish_mode)

        combineClinvarGenesResolvedCalls(combineResolvedCalls.out, lowGenes)
        countClinvarGeneResolvedCalls(combineClinvarGenesResolvedCalls.out, workflow.projectDir)
        addParamsToClinvarResolvedCallsCombined(countClinvarGeneResolvedCalls.out, filterStr, workflow.projectDir, publish_mode)
}
