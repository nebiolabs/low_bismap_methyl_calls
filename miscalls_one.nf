nextflow.enable.dsl=2

params.clinvarBed="test_fixtures/chr18_clinvar_regions.bed"

params.bismap="test_fixtures/chr18_bismap.bw"

params.bismapBbm="test_fixtures/chr18_bismap.bbm"

params.refFai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams=""

params.tmpdir = "/tmp"

params.minMapq=10

params.noFilterMappability=false

params.bismapCutoff=0.01

params.hardlink = false

params.outputDir="methyl_calls"

params.trimBases=8

params.pureCutoff=10

outputDir=params.outputDir

process bigWigToBedGraph {

        penv 'smp'
        cpus 1
        conda 'ucsc-bigwigtobedgraph=377'

        input:
                file bismap_bw

        output:
                file 'bismap.bedGraph'

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bigWigToBedGraph !{bismap_bw} bismap.bedGraph
                '''
}

process convertFaiToGenome {

        penv 'smp'
        cpus 1
        conda 'coreutils=9.1'
        input:
                file ref_fai
                val tmpdir

        output:
                file 'ref.genome'

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                cut -f 1,2 !{ref_fai} | LC_ALL=C sort -k1,1 -T !{tmpdir}  > ref.genome
                '''
}

process extractBam {

        penv 'smp'
        conda 'methyldackel zstd=1.5.6'
        cpus 8

        input:
                file ref
                val min_mapq
                tuple val(_groupKey), file(bam), file(bai)
                val filter_opt
                file ref_fai
                file bismap_bbm
  
        output:
                tuple file('*_CpG.bedGraph.zst'), file('*_CHG.bedGraph.zst'), file('*_CHH.bedGraph.zst'), val("${bam.getBaseName().replace(".bam\n", "")}")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                MethylDackel extract -@ !{task.cpus} -q !{min_mapq} --CHH --CHG --nOT 0,0,0,!{params.trimBases} --nOB 0,0,!{params.trimBases},0 !{filter_opt} -o methyl_!{bam.getBaseName().replace(".bam\\n", "")} !{ref} !{bam}
                for bedGraph in *.bedGraph; do
                    context="$(basename "$bedGraph" | sed -E "s/methyl_!{bam.getBaseName().replace(".bam\\n", "")}_([CGHp]+)\\.bedGraph/\\1/")"
                    sed -E "s/$/\t$context/" -i $bedGraph;
                    zstd --rm "$bedGraph" -o "$bedGraph.zst";
                done
                '''
}

process extractBamMq0 {

        penv 'smp'
        conda 'methyldackel zstd=1.5.6'
        cpus 8

        input:
                file ref
                val min_mapq
                tuple val(_groupKey), file(bam), file(bai)
                val filter_opt
                file ref_fai
                file bismap_bbm
  
        output:
                tuple file('*_CpG.bedGraph.zst'), file('*_CHG.bedGraph.zst'), file('*_CHH.bedGraph.zst'), val("mq0_${bam.getBaseName().replace(".bam\n", "")}")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                MethylDackel extract -@ !{task.cpus} -q 0 --CHH --CHG --nOT 0,0,0,!{params.trimBases} --nOB 0,0,!{params.trimBases},0 !{filter_opt} -o methyl_mq0_!{bam.getBaseName().replace(".bam\\n", "")} !{ref} !{bam}
                
                for bedGraph in *.bedGraph; do
                    context="$(basename "$bedGraph" | sed -E "s/methyl_mq0_!{bam.getBaseName().replace(".bam\\n", "")}_([CGHp]+)\\.bedGraph/\\1/")"
                    sed -E "s/$/\t$context/" -i $bedGraph;
                    zstd --rm "$bedGraph" -o "$bedGraph.zst";
                done
                '''
}
process combineCalls {

        conda 'coreutils=9.1 zstd=1.5.6'
        penv 'smp'
        cpus 8

        input:
                tuple file(methyl_cpg), file(methyl_chg), file(methyl_chh), val(n)
                val tmpdir

        output:
                tuple file('sorted_noheader_*.bedGraph.zst'), val(n)
    
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} <(zstd -d -f --stdout !{methyl_cpg}) <(zstd -d -f --stdout !{methyl_chg}) <(zstd -d -f --stdout !{methyl_chh}) | head -n -3 | zstd -o sorted_noheader_!{n.strip()}.bedGraph.zst
                '''
}

process complementBedgraph {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0'

        input:
                file ref_genome
                file bismap_bedGraph_for_complement

        output:
                file 'bismap_zeroes.bedGraph'
    
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools complement -i !{bismap_bedGraph_for_complement} -g !{ref_genome} | awk '{ print $0"\t"0 }'> bismap_zeroes.bedGraph
                '''
}

process combineAllRegions {

        penv 'smp'
        cpus 8
        conda 'coreutils=9.1'

        input:
                file bismap_bedGraph_for_combine_all
                file zeroes_bedgraph
                val tmpdir

        output:
                file 'bismap_all.bedGraph'
    
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{bismap_bedGraph_for_combine_all} !{zeroes_bedgraph} > bismap_all.bedGraph
                '''
}

process intersectMin {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0'

        input:
                file clinvar_regions_for_min
                file bismap_bedgraph_for_min

        output:
                file 'clinvar_regions_min.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools map -o min -c 4 -a !{clinvar_regions_for_min} -b !{bismap_bedgraph_for_min} > clinvar_regions_min.bedGraph
                '''
}

process intersectMean {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0'

        input:
                file clinvar_regions_for_mean
                file bismap_bedgraph_for_mean

        output:
                file 'clinvar_regions_mean.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools map -o mean -c 4 -a !{clinvar_regions_for_mean} -b !{bismap_bedgraph_for_mean} > clinvar_regions_mean.bedGraph
                '''
}

process pasteMinMean {

        penv 'smp'
        cpus 1
        conda 'coreutils=9.1'
        input:
                file clinvar_regions_mean_map
                file clinvar_regions_min_map

        output:
                file 'clinvar_regions_min_mean.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                paste !{clinvar_regions_min_map} <(cut -f 11 !{clinvar_regions_mean_map}) > clinvar_regions_min_mean.bedGraph
                '''
}

process filterMinMean {

        penv 'smp'
        cpus 1
        conda 'gawk=5.3.0'
        input:
                file clinvar_regions_min_mean_map
                val bismap_cutoff

        output:
                file 'clinvar_regions_min_mean_low.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                awk '$11 <= !{bismap_cutoff} { print }' !{clinvar_regions_min_mean_map} > clinvar_regions_min_mean_low.bedGraph
                '''
}

process gsortMinMean {

        cpus 8
        penv 'smp'
        conda 'coreutils=9.1'

        input:
                file clinvar_regions_min_mean_map_low_unsorted
                val tmpdir

        output:
                file 'clinvar_annotations_low.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{clinvar_regions_min_mean_map_low_unsorted} > clinvar_annotations_low.bedGraph
                '''
}

process slopRegions {

        cpus 8
        penv 'smp'
        conda 'bedtools=2.30.0 coreutils=9.1 gawk=5.3.0'

        input:
                file clinvar_regions_min_mean_map_low
                val tmpdir
                file ref_genome

        output:
                file 'clinvar_annotations_low_expanded_unsorted.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                cat !{clinvar_regions_min_mean_map_low} | awk '{ print $0"\tplaceholder" }' > tmp.bed; bedtools slop -l 2000 -r 0 -s -i tmp.bed -g !{ref_genome} | cut -f 1-12 > clinvar_annotations_low_expanded_unsorted.bedGraph; rm tmp.bed
                '''
}

process gsortSlopRegions {

        cpus 8
        penv 'smp'
        conda 'coreutils=9.1'

        input:
                file slop_regions
                val tmpdir

        output:
                file 'clinvar_annotations_low_expanded.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{slop_regions} > clinvar_annotations_low_expanded.bedGraph
                '''
}


process getLowRegions {

        penv 'smp'
        cpus 1
        conda 'gawk=5.3.0'
        input:
                file bismap_bedGraph_for_awk_low_regions

        output:
                file 'bismap_low_nozeroes.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                awk '$4 <= 0.01 { print }' !{bismap_bedGraph_for_awk_low_regions} > bismap_low_nozeroes.bedGraph
                '''
}



process combineLowRegions {

        cpus 8
        penv 'smp'
        conda 'coreutils=9.1'

        input:
                file low_bedgraph
                file zeroes_bedgraph
                val tmpdir

        output:
                file 'bismap_low.bedGraph'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                LC_ALL=C sort -k1,1 -k2,2n -S 32G -T !{tmpdir} --parallel !{task.cpus} !{low_bedgraph} !{zeroes_bedgraph} > bismap_low.bedGraph
                '''
}

process intersectClinvarBismapLow {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0'

        input:
                file clinvar_annotations_low
                file bismap_low

        output:
                file 'clinvar_low_regions.bed'
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -sorted -a !{bismap_low} -b !{clinvar_annotations_low} > clinvar_low_regions.bed
                '''
}

process intersectCalls {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 gawk=5.3.0 zstd=1.5.6'

        input:
                file regions_for_methyl
                tuple file(methyl_sorted), val(n)

        output:
                tuple file("lowmap_calls_*.bed.zst"), val(n)

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -u -sorted -a <(zstd -d -f --stdout !{methyl_sorted}) -b !{regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' | zstd -o lowmap_calls_!{n.strip()}.bed.zst
                '''
}

process annotateCalls {

        penv 'smp'
        cpus 1
        conda 'coreutils=9.1 gawk=5.3.0 zstd=1.5.6'
        input:
                tuple file(methyl_sorted), val(n)

        output:
                tuple file("annotated_calls_*.bed.zst"), val(n)

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{methyl_sorted} | awk '{ print $0"\t!{n.strip()}" }' | zstd -o annotated_calls_!{n.strip()}.bed.zst
                '''
}

process intersectClinvarCalls {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 gawk=5.3.0 zstd=1.5.6'

        input:
                file clinvar_regions_for_methyl
                tuple file(methyl_sorted), val(n)

        output:
                tuple file("clinvar_calls_*.bed.zst"), val(n)

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -u -sorted -a <(zstd -d -f --stdout !{methyl_sorted}) -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' | zstd -o clinvar_calls_!{n.strip()}.bed.zst
                '''
}

process intersectClinvarCallsLowmap {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 gawk=5.3.0 zstd=1.5.6'

        input:
                file clinvar_regions_for_methyl
                tuple file(methyl_sorted), val(n)

        output:
                tuple file("clinvar_lowmap_calls_*.bed.zst"), val(n)

        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -u -sorted -a <(zstd -d -f --stdout !{methyl_sorted}) -b !{clinvar_regions_for_methyl} | awk '{ print $0"\t!{n.strip()}" }' | zstd -o clinvar_lowmap_calls_!{n.strip()}.bed.zst
                '''
}


process intersectCalledGenes {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple file(methyl_calls), val(n)
                file clinvar_annotations_low

        output:
                tuple file('called_genes_dups_*.bed.zst'), val(n)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -wa -wb -sorted -a !{clinvar_annotations_low} -b <(zstd -d -f --stdout !{methyl_calls}) | zstd -o called_genes_dups_!{n.strip()}.bed.zst
                '''
}


process intersectLowmapCalledGenes {

        penv 'smp'
        cpus 1
        conda 'bedtools=2.30.0 zstd=1.5.6'

        input:
                tuple file(methyl_lowmap_calls), val(n)
                file clinvar_annotations_low 

        output:
                tuple file('lowmap_called_genes_dups_*.bed.zst'), val(n)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                bedtools intersect -wa -wb -sorted -a !{clinvar_annotations_low} -b <(zstd -d -f --stdout !{methyl_lowmap_calls}) | zstd -o lowmap_called_genes_dups_!{n.strip()}.bed.zst
                '''
}

process countLowmapGeneCalls {

        penv 'smp'
        cpus 1
        conda 'python=3.12.2 zstd=1.5.6'
        //use nextflow bin dir for script?
        input:
                tuple file(lowmap_called_genes_dups_methyl), val(n)
                val basePath

        output:
                tuple file("lowmap_called_genes_*.bed.zst"), val(n)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_called_genes_dups_methyl} | python !{basePath}/count_calls.py !{n.strip()} | zstd -o lowmap_called_genes_!{n.strip()}.bed.zst 
                '''
}

process countGeneCalls {

        penv 'smp'
        cpus 1
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(called_genes_dups_methyl), val(n)
                val basePath

        output:
                tuple file("called_genes_*.bed.zst"), val(n)
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{called_genes_dups_methyl} | python !{basePath}/count_calls.py !{n.strip()} | zstd -o called_genes_!{n.strip()}.bed.zst
                '''
}


process addParamsToLowmapCalledGenes {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/lowmap_called_genes", mode: "${publish_mode}", pattern: "lowmap_called_genes_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(lowmap_called_genes), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("lowmap_called_genes_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_called_genes} | python !{basePath}/add_params.py !{filterStr} 15 | zstd -o lowmap_called_genes_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process addParamsToCalledGenes {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/called_genes", mode: "${publish_mode}", pattern: "called_genes_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(called_genes), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("called_genes_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{called_genes} | python !{basePath}/add_params.py !{filterStr} 15 | zstd -o called_genes_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process addParamsToCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/calls", mode: "${publish_mode}", pattern: "calls_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(calls), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("calls_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{calls} | python !{basePath}/add_params.py !{filterStr} 7 | zstd -o calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process addParamsToLowmapCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/lowmap_calls", mode: "${publish_mode}", pattern: "lowmap_calls_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(lowmap_calls), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("lowmap_calls_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_calls} | python !{basePath}/add_params.py !{filterStr} 7 | zstd -o lowmap_calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process addParamsToClinvarCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/clinvar_calls", mode: "${publish_mode}", pattern: "clinvar_calls_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(lowmap_calls), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("clinvar_calls_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_calls} | python !{basePath}/add_params.py !{filterStr} 7 | zstd -o clinvar_calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process addParamsToClinvarLowmapCalls {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/clinvar_lowmap_calls", mode: "${publish_mode}", pattern: "clinvar_lowmap_calls_*.${filterStr}.tsv.zst"
        conda 'python=3.12.2 zstd=1.5.6'

        input:
                tuple file(lowmap_calls), val(n)
                val filterStr
                val basePath
                val publish_mode

        output:
                tuple val(n), file("clinvar_lowmap_calls_*.${filterStr}.tsv.zst")
        
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                zstd -d -f --stdout !{lowmap_calls} | python !{basePath}/add_params.py !{filterStr} 7 | zstd -o clinvar_lowmap_calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process getPureMethylSites {

        penv 'smp'
        cpus 1
        publishDir "${outputDir}/pure_mixed_calls", mode: "${publish_mode}", pattern: "pure_methyl_calls_*.${filterStr}.tsv.zst"
        conda 'gawk=5.3.0 zstd=1.5.6'

        input:
                tuple val(n), file(calls)
                val filterStr
                val publish_mode
        output:
                tuple val(n), file("pure_methyl_calls_*.${filterStr}.tsv.zst")
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                awk '($4 >= (100-!{params.pureCutoff}) || $4 <= !{params.pureCutoff}) { print }' <(zstd -d -f --stdout !{calls}) | zstd -o pure_methyl_calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

process getMixedMethylSites {

        cpus 2
        penv 'smp'
        publishDir "${outputDir}/pure_mixed_calls", mode: "${publish_mode}", pattern: "mixed_methyl_calls_*.${filterStr}.tsv.zst"
        conda 'gawk=5.3.0 zstd=1.5.6'

        input:
                tuple val(n), file(calls)
                val filterStr
                val publish_mode
        output:
                tuple val(n), file("mixed_methyl_calls_*.${filterStr}.tsv.zst")
        shell:
                '''
                export TMPDIR="!{params.tmpdir}"
                awk '($4 < (100-!{params.pureCutoff}) && $4 > !{params.pureCutoff}) { print }' <(zstd -d -f --stdout !{calls}) | zstd -o mixed_methyl_calls_!{n.strip()}.!{filterStr}.tsv.zst
                '''
}

workflow {

    no_filter_mappability=params.noFilterMappability

    bismap_bw = file(params.bismap)
    clinvar_regions = file(params.clinvarBed)
    ref_fai = file(params.refFai)
    ref = file(params.ref)
    bismap_bbm= file(params.bismapBbm)

    pipeline_bams = Channel.fromFilePairs(params.bam_glob+"{.bam,.bam.bai}")

    publishMode = "symlink"

    if(params.hardlink)
    {
        publishMode = "link"
    }

    refPreproc(bismap_bw, ref, ref_fai, clinvar_regions, params.tmpdir, params.minMapq, params.bismapCutoff)
    
    basicPipeline(bismap_bbm, ref, ref_fai, no_filter_mappability, pipeline_bams, params.tmpdir, params.minMapq, params.bismapCutoff, refPreproc.out[0], refPreproc.out[1], refPreproc.out[2], publishMode)


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
        slopRegions(combineAllRegions.out, tmpdir, convertFaiToGenome.out)
        gsortSlopRegions(slopRegions.out, tmpdir)
        getLowRegions(bigWigToBedGraph.out)
        combineLowRegions(getLowRegions.out, complementBedgraph.out, tmpdir)
        intersectMin(clinvar_regions, gsortSlopRegions.out)
        intersectMean(clinvar_regions, gsortSlopRegions.out)
        pasteMinMean(intersectMean.out, intersectMin.out)
        filterMinMean(pasteMinMean.out, bismap_cutoff)
        gsortMinMean(filterMinMean.out, tmpdir)
        intersectClinvarBismapLow(gsortMinMean.out, combineLowRegions.out)
    emit:
        intersectClinvarBismapLow.out
        combineLowRegions.out
        gsortMinMean.out
}

workflow basicPipeline {
    take:
        bismap_bbm
        ref
        ref_fai
        no_filter_mappability
        pipeline_bams
        tmpdir
        min_mapq
        bismap_cutoff
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

        //output Suffix = referenceName+"_"+aligner+"_"+seqMethod+"_"+filterStr

        //bam_names = pipeline_bams.toString().split(",")

        //println(bam_names)

        //bam_paths = bam_names.collect({ "readlink -f ${it}".execute().text })
        //bam_list = Channel.from(bam_paths)

        extractBam(ref, min_mapq, pipeline_bams, filter_opt, ref_fai, bismap_bbm)
        extractBamMq0(ref, min_mapq, pipeline_bams, filter_opt, ref_fai, bismap_bbm)
        combineCalls(extractBam.out.mix(extractBamMq0.out), tmpdir)
        intersectClinvarCalls(lowGenes, combineCalls.out)
        intersectClinvarCallsLowmap(clinvarBismapLow, combineCalls.out)
        intersectCalls(lowRegions, combineCalls.out)
        annotateCalls(combineCalls.out)
        intersectLowmapCalledGenes(intersectCalls.out, lowGenes)
        intersectCalledGenes(annotateCalls.out, lowGenes)
        countLowmapGeneCalls(intersectLowmapCalledGenes.out, workflow.projectDir)
        countGeneCalls(intersectCalledGenes.out, workflow.projectDir)
        addParamsToCalls(annotateCalls.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToLowmapCalls(intersectCalls.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToClinvarCalls(intersectClinvarCalls.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToClinvarLowmapCalls(intersectClinvarCallsLowmap.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToCalledGenes(countGeneCalls.out, filterStr, workflow.projectDir, publish_mode)
        addParamsToLowmapCalledGenes(countLowmapGeneCalls.out, filterStr, workflow.projectDir, publish_mode)
        getPureMethylSites(addParamsToCalls.out, filterStr, publish_mode)
        getMixedMethylSites(addParamsToCalls.out, filterStr, publish_mode)
    emit:
        addParamsToLowmapCalls.out
        addParamsToCalls.out
        addParamsToClinvarLowmapCalls.out
        getPureMethylSites.out
        getMixedMethylSites.out
        addParamsToClinvarCalls.out
}
