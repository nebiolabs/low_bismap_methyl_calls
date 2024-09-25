nextflow.enable.dsl=2

include { get_unique_calls as get_unique_calls_emseq; get_unique_calls as get_unique_miscalls_emseq; get_unique_calls as get_unique_calls_wgbs; get_unique_calls as get_unique_miscalls_wgbs } from './get_unique_calls.nf'

process get_lib_params {
    penv 'smp'
    cpus 1
   
    input:
        val dir
        val tmp_dir
 
    output:
        file "*.txt"

    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    for lib_params in $(find !{dir} -mindepth 1 -maxdepth 1 -type d); do
        echo "$lib_params" > "$(basename "${lib_params}")".txt
    done 
    '''
}

process calc_common_calls {
    penv 'smp'
    cpus 1
    publishDir "${output_dir}/${lib_params.getBaseName()}", mode: "${publish_mode}", pattern: "*.bed.gz"
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20 zstd=1.5.6'    
   
    input:
        file lib_params
        val output_dir
        val publish_mode
        val tmp_dir
 
    output:
        tuple file("*.bed.gz"), val("${lib_params.getBaseName()}")
    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    dir_path="$(cat "!{lib_params}")"
    bedtools intersect -sorted -a <(zstdcat ${dir_path}/emseq/calls_nofilt.bed.zst) -b <(zstdcat ${dir_path}/wgbs/calls_nofilt.bed.zst) | bgzip -o "common.calls.!{lib_params.getBaseName()}.bed.gz"
    '''
}

process calc_common_miscalls {
    penv 'smp'
    cpus 1
    publishDir "${output_dir}/${lib_params.getBaseName()}", mode: "${publish_mode}", pattern: "*.bed.gz"
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20 zstd=1.5.6'    
   
    input:
        file lib_params
        val output_dir
        val publish_mode
        val tmp_dir
 
    output:
        tuple file("*.bed.gz"), val("${lib_params.getBaseName()}")
    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    dir_path="$(cat "!{lib_params}")"
    bedtools intersect -sorted -u -a <(zstdcat ${dir_path}/emseq/miscalls.bed.zst) -b <(zstdcat ${dir_path}/wgbs/miscalls.bed.zst) | bgzip -o "common.miscalls.!{lib_params.getBaseName()}.bed.gz"
    '''
}

process ensure_completed {
    penv 'smp'
    cpus 1

    input:
        val done
        val output_dir

    output:
        env outdir

    shell:
    '''
    export outdir="!{output_dir}"
    '''
}

workflow processRef {

    take:
        libs_dir
        output_dir
        publish_mode
        tmp_dir

    main:
        get_lib_params(libs_dir, tmp_dir)
        lib_params = get_lib_params.out.flatten()
        calc_common_calls(lib_params, output_dir, publish_mode, tmp_dir)
        calc_common_miscalls(lib_params, output_dir, publish_mode, tmp_dir)
        get_unique_calls_emseq(lib_params, "emseq", "wgbs", "calls_nofilt", output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_emseq(lib_params, "emseq", "wgbs", "miscalls", output_dir, publish_mode, tmp_dir)
        get_unique_calls_wgbs(lib_params, "wgbs", "emseq", "calls_nofilt", output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_wgbs(lib_params, "wgbs", "emseq", "miscalls", output_dir, publish_mode, tmp_dir)
        unique_calls = get_unique_calls_emseq.out.mix(get_unique_miscalls_emseq.out, get_unique_calls_wgbs.out, get_unique_miscalls_wgbs.out)
        ensure_completed(unique_calls.collect(), output_dir) 

    emit:
        ensure_completed.out

}

workflow {
    processRef(params.libsDir, params.bismapBw, params.chromSizes, params.outputDir, params.publishMode, params.tmpdir)
}

