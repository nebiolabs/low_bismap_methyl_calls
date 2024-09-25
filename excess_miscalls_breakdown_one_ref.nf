nextflow.enable.dsl=2

include { get_unique_miscalls as get_unique_miscalls_emseq; get_unique_miscalls as get_unique_miscalls_emseq_clinvar; get_unique_miscalls as get_unique_miscalls_wgbs; get_unique_miscalls as get_unique_miscalls_wgbs_clinvar } from './get_unique_miscalls.nf'

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
    bedtools intersect -sorted -u -a <(zstdcat ${dir_path}/emseq/clinvar_miscalls.bed.zst) -b <(zstdcat ${dir_path}/wgbs/clinvar_miscalls.bed.zst) | bgzip -o "common.clinvar_miscalls.!{lib_params.getBaseName()}.bed.gz"
    '''
}

process get_bismap_bedgraph {
    penv 'smp'
    cpus 1
    conda 'ucsc-bigwigtobedgraph=377'    

    input:
        file bismap_bigwig
        val tmp_dir
 
    output:
        file "*bedGraph"

    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    bigWigToBedGraph "!{bismap_bigwig}" bismap_bedgraph.bedGraph
    '''
}

process get_common_unique_miscalls {
    penv 'smp'
    cpus 1
    publishDir "${output_dir}/${lib_params.getBaseName()}", mode: "${publish_mode}", pattern: "*.bed.gz"
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20'
    
    input:
        tuple val(_joinKey), file(unique_miscalls), file(lib_params), file(common_calls)
        val output_dir
        val publish_mode
        val tmp_dir
    
    output:
        tuple file("varying.shared.*.bed.gz"), file(lib_params)
        file "varying.unique.*.bed.gz"

    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    bedtools intersect -sorted -wa -u -a <(zcat !{unique_miscalls}) -b <(zcat !{common_calls}) | bgzip -o varying.shared.!{unique_miscalls}
    bedtools intersect -v -sorted -wa -a <(zcat !{unique_miscalls}) -b <(zcat !{common_calls}) | bgzip -o varying.unique.!{unique_miscalls}
    '''
}

process get_near_mappable_regions {
    penv 'smp'
    cpus 1
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20'    

    input:
        file bismap_bedGraph
        file chrom_sizes
        val tmp_dir
 
    output:
        file "*bed.gz"

    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    bedtools merge -i !{bismap_bedGraph} | bedtools slop -i stdin -g !{chrom_sizes} -b 500 | bgzip -o near_mappable_regions.bed.gz
    '''
}

process intersect_near_mappable_regions {
    penv 'smp'
    cpus 1
    publishDir "${output_dir}/${lib_params.getBaseName()}", mode: "${publish_mode}", pattern: "intersected.*.bed.gz"
    conda 'bedtools=2.30.0 tabix=1.11 htslib=1.20'    
   
    input:
        tuple file(extra_miscalls), file(lib_params), file(near_mappable_regions)
        val output_dir
        val publish_mode
        val tmp_dir

    output:
        file "intersected.*.bed.gz"
        val "done"

    shell:
    '''
    export TMPDIR="!{tmp_dir}"
    bedtools intersect -sorted -wa -u -a <(zcat !{extra_miscalls}) -b <(zcat !{near_mappable_regions}) | bgzip -o "intersected.close.!{extra_miscalls}"
    bedtools intersect -sorted -v -wa -a <(zcat !{extra_miscalls}) -b <(zcat !{near_mappable_regions}) | bgzip -o "intersected.far.!{extra_miscalls}"
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
        bismap_bw
        chrom_sizes
        output_dir
        publish_mode
        tmp_dir

    main:
        get_lib_params(libs_dir, tmp_dir)
        get_bismap_bedgraph(file(bismap_bw), tmp_dir)
        get_near_mappable_regions(get_bismap_bedgraph.out, file(chrom_sizes), tmp_dir)
        lib_params = get_lib_params.out.flatten()
        calc_common_calls(lib_params, output_dir, publish_mode, tmp_dir)
        calc_common_miscalls(lib_params, output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_emseq(lib_params, "emseq", "wgbs", "miscalls", output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_emseq_clinvar(lib_params, "emseq", "wgbs", "clinvar_miscalls", output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_wgbs(lib_params, "wgbs", "emseq", "miscalls", output_dir, publish_mode, tmp_dir)
        get_unique_miscalls_wgbs_clinvar(lib_params, "wgbs", "emseq", "clinvar_miscalls", output_dir, publish_mode, tmp_dir)
        unique_miscalls = get_unique_miscalls_emseq.out.mix(get_unique_miscalls_emseq_clinvar.out, get_unique_miscalls_wgbs.out, get_unique_miscalls_wgbs_clinvar.out)
        get_common_unique_miscalls(unique_miscalls.combine(calc_common_calls.out, by: 1), output_dir, publish_mode, tmp_dir)
        intersect_near_mappable_regions(get_common_unique_miscalls.out[0].combine(get_near_mappable_regions.out), output_dir, publish_mode, tmp_dir)
        ensure_completed(intersect_near_mappable_regions.out[1].collect(), output_dir) 

    emit:
        ensure_completed.out

}

workflow {
    processRef(params.libsDir, params.bismapBw, params.chromSizes, params.outputDir, params.publishMode, params.tmpdir)
}

