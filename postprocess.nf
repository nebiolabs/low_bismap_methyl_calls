process count_per_chrom {
    penv 'smp'
    cpus 1
    publishDir "${params.outputDir}/counts", mode: "${params.publishMode}", pattern: "*.tsv"
    conda 'zstd=1.5.6'    
   
    input:
        file data_tsv    
 
    output:
        file "*.tsv"

    shell:
    '''
    export TMPDIR="!{params.tmpdir}"
    zstdcat !{data_tsv} | cut -f 1 | uniq -c | sed -E 's/^\\s+//; s/ /\t/' > "counts_!{data_tsv.getBaseName()}"
    '''
}

process summarize_counts {
    penv 'smp'
    cpus 1
    publishDir "${params.outputDir}/counts_summary", mode: "${params.publishMode}", pattern: "totals_*.tsv"    
   
    input:
        file counts_tsvs   
 
    output:
        file "totals_*.tsv"

    shell:
    '''
    export TMPDIR="!{params.tmpdir}"
    while read -d " " -r tsv; do
        prefixes="called_genes|calls|clinvar_calls|clinvar_lowmap_calls|clinvar_miscalls|clinvar_miscalls_lowmap|lowmap_called_genes|lowmap_calls|miscalled_clinvar_genes|miscalled_clinvar_genes_lowmap|miscalls|miscalls_lowmap|mixed_methyl_calls|pure_methyl_calls|resolved_calls|resolved_clinvar_genes|underfiltered_calls"
        metadata="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_([S0-9][^.]+)\\.([^.]+)\\.(bisulfite|emseq)\\.(grch38|t2t)\\.(bismark|bwameth)(\\.(filt|nofilt))?\\.tsv/\\4.\\5.\\6.\\7.\\8/")"
        processing_stage="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_.*.tsv/\\1/")"
        filt_status="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_[S0-9][^.]+\\.[^.]+\\.[^.]+\\.[^.]+\\.[^.]+(\\.(filt|nofilt))?\\.tsv/\\5/")"
        output_metadata="$(echo "$tsv" | sed -E "s/counts_($prefixes)(_mq0)?_([S0-9][^.]+)\\.([^.]+)\\.(bisulfite|emseq)\\.(grch38|t2t)\\.(bismark|bwameth).*\\.tsv/\\3\t\\4\t\\5\t\\6\t\\7/")"
        if [ -z "$filt_status" ]; then 
            filt_status="combined";
        fi
        output_metadata="$filt_status\t$output_metadata"
        total="$(cut -f 1 "$tsv" | paste -sd+ | bc)" 
        echo "$total\t${processing_stage}\t${output_metadata}" >> "totals_${metadata}.tsv"
    done <<< "!{counts_tsvs} "
    '''
}

process sankey {
    penv 'smp'
    cpus 1
    conda "python plotly python-kaleido"
    publishDir "${params.outputDir}", mode: "${params.publishMode}", pattern: "sankey_*.svg"
    
    input:
        file lib_summary

    output:
        file "sankey_*.svg"

    shell:
    '''
    export TMPDIR="!{params.tmpdir}"
    python !{workflow.projectDir}/plotly_sankey.py "!{lib_summary}" "!{lib_summary.getBaseName().replace("totals_", "").replace(".tsv", "")}"
    '''
}

workflow {

    calls = Channel.fromPath(params.methylCallsDir+"/calls/*.tsv.zst")
    clinvar_calls = Channel.fromPath(params.methylCallsDir+"/clinvar_calls/*.tsv.zst")
    clinvar_lowmap_calls = Channel.fromPath(params.methylCallsDir+"/clinvar_lowmap_calls/*.tsv.zst")
    lowmap_calls = Channel.fromPath(params.methylCallsDir+"/lowmap_calls/*.tsv.zst")
    miscalls = Channel.fromPath(params.methylCallsDir+"/miscalls/*.tsv.zst")
    pure_mixed_calls = Channel.fromPath(params.methylCallsDir+"/pure_mixed_calls/*.tsv.zst")
    underfiltered_calls = Channel.fromPath(params.methylCallsDir+"/underfiltered_calls/*.tsv.zst")
    called_genes = Channel.fromPath(params.methylCallsDir+"/called_genes/*.tsv.zst")
    lowmap_called_genes = Channel.fromPath(params.methylCallsDir+"/lowmap_called_genes/*.tsv.zst")
    miscalled_genes = Channel.fromPath(params.methylCallsDir+"/miscalled_genes/*.tsv.zst")
    
    all_calls = calls.mix(clinvar_calls, clinvar_lowmap_calls, lowmap_calls, miscalls, pure_mixed_calls, underfiltered_calls, called_genes, lowmap_called_genes, miscalled_genes)

    count_per_chrom(all_calls)
    summarize_counts(count_per_chrom.out.collect())
    sankey(summarize_counts.out.flatten())
}

