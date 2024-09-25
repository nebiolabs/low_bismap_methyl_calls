nextflow.enable.dsl=2

process get_unique_miscalls {
    penv 'smp'
    cpus 1
    publishDir "${output_dir}/${lib_params.getBaseName()}", mode: "${publish_mode}", pattern: "*.bed.gz"
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20 zstd=1.5.6'

    input:
        file lib_params
        val method_1
        val method_2
        val bed_name
        val output_dir
        val publish_mode
        val tmpdir

    output:
        tuple file("*bed.gz"), val("${lib_params.getBaseName()}"), file(lib_params)

    shell:
    '''
    dir_path="$(cat "!{lib_params}")"
    export TMPDIR="!{tmpdir}"
    bedtools subtract -sorted -a <(zstdcat ${dir_path}/!{method_1}/!{bed_name}.bed.zst) -b <(zstdcat ${dir_path}/!{method_2}/!{bed_name}.bed.zst) | awk 'BEGIN { OFS="\\t" } { print $1, $2, $3, $5+$6, $4, $5, $6, $7, $8, $9, $10, $11 }'  | bgzip -o excess.!{method_1}.!{bed_name}.bed.gz
    '''
}
