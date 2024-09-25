nextflow.enable.dsl=2

include { processRef as processRefGRCh38 } from './calls_breakdown_single.nf'
include { processRef as processRefT2T } from './calls_breakdown_single.nf'

process countCategorizedMiscalls {
    penv 'smp'
    cpus 1
    publishDir "${params.outputDir}", mode: "${params.publishMode}", pattern: "counts_summary.tsv"
    conda 'gzip=1.13'

    input:
        val data_dir_grch38
        val data_dir_t2t
    
    output:
        file "counts_summary.tsv"

    shell:
    '''
    export TMPDIR="!{params.tmpdir}"
    for outdir in !{data_dir_grch38} !{data_dir_t2t}; do
        for dir in $outdir/*; do
            dirinfo="$(basename "$dir" | tr "-" "\t")";
            for bed in $dir/common.*.bed.gz; do
                bedname="$(basename "$bed")"
                echo "$(zcat "$bed" | wc -l)\t$(echo "$bedname" | cut -f 1-2 -d ".")\t$dirinfo" >> "counts_summary.tsv"
            done
            for bed in $dir/excess.*.bed.gz; do
                bedname="$(basename "$bed")"
                echo "$(zcat "$bed" | wc -l)\t$(echo "$bedname" | cut -f 1-3 -d ".")\t$dirinfo" >> "counts_summary.tsv"
            done
        done
    done
    '''
}


workflow {
    processRefGRCh38(params.libsDirGrch38, params.outputDirGrch38, params.publishMode, params.tmpdir)
    processRefT2T(params.libsDirT2t, params.outputDirT2t, params.publishMode, params.tmpdir)
    countCategorizedMiscalls(processRefGRCh38.out, processRefT2T.out)
}
