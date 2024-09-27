nextflow.enable.dsl=2

assert params.output_dir != null : "--output_dir must be specified"
assert params.tmp_dir != null : "--tmp_dir must be specified"

process add_rg {
    cpus 1

    conda "samtools=1.19"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    file input_bam

    output:
    file "*.rg.bam"

    shell:
    '''
    export TMPDIR="!{params.tmp_dir}"
    samtools addreplacerg -r "@RG\tID:-\tSM:-" -o "$(basename -s .bam "!{input_bam}").rg.bam" "!{input_bam}" 
    '''
    }

process index_bams {
    cpus 1

    conda "samtools=1.19"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    file input_bam

    output:
    file "*.bai"

    shell:
    '''
    export TMPDIR="!{params.tmp_dir}"
    samtools index -o "!{input_bam}.bai" "!{input_bam}"
    '''
    }

workflow {
    bams = Channel.fromPath( '*.bam' )
    add_rg(bams)
    index_bams(add_rg.out)
}
