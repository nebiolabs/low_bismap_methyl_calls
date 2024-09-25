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
    zstdcat !{data_tsv} | awk 'BEGIN { FS="\t"; OFS="\t"; } /\tCpG/ { count_cpg[$1]+=1; next; } /\tCHG/ { count_chg[$1]+=1; next; } /\tCHH/ { count_chh[$1]+=1; next; } END { for(chr in count_cpg) { print count_cpg[chr], chr, "CpG" }; for(chr in count_chg) { print count_chg[chr], chr, "CHG" }; for(chr in count_chh) { print count_chh[chr], chr, "CHH" }; }' > "counts_!{data_tsv.getBaseName()}"
    '''
}

process count_per_chrom_genes {
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
    zstdcat !{data_tsv} | awk 'BEGIN { FS="\t"; OFS="\t"; } $13 > 0 { count_cpg[$1]+=1; } $14 > 0 { count_chg[$1]+=1; } $15 > 0 { count_chh[$1]+=1; next; } END { for(chr in count_cpg) { print count_cpg[chr], chr, "CpG" }; for(chr in count_chg) { print count_chg[chr], chr, "CHG" }; for(chr in count_chh) { print count_chh[chr], chr, "CHH" }; }' > "counts_!{data_tsv.getBaseName()}"
    '''
}

process count_ref_cs {
    penv 'smp'
    cpus 1
    conda 'coreutils seqtk'
    
    input:
        file seq

    output:
        file "seq_counts_${seq.getBaseName()}.tsv"

    shell:
    '''
    count=$(grep -v ">" !{seq} | tr -c -d 'cCgG' | wc -c)
    count_cpg_1strand=$(seqtk seq !{seq} | grep -v ">" | grep -o -E -i "cg" | wc -l)
    count_cpg=$((2* count_cpg_1strand))
    echo "!{seq.getBaseName().replace(".fasta", "").replace(".fa", "")}_all\t$count" >> "seq_counts_!{seq.getBaseName()}.tsv"
    echo "!{seq.getBaseName().replace(".fasta", "").replace(".fa", "")}_cpg\t$count_cpg" >> "seq_counts_!{seq.getBaseName()}.tsv"
    '''
}

process combine_c_counts {
    penv 'smp'
    cpus 1
    
    input:
        file counts_tsvs

    output:
        file "c_counts.tsv"

    shell:
    '''
    cat !{counts_tsvs} > c_counts.tsv
    '''
}

process summarize_counts {
    penv 'smp'
    cpus 1
    publishDir "${params.outputDir}/counts_summary", mode: "${params.publishMode}", pattern: "totals_*.tsv"    
   
    input:
        file counts_tsvs 
        file c_counts
 
    output:
        file "totals_*.tsv"

    shell:
    '''
    export TMPDIR="!{params.tmpdir}"
    while read -d " " -r tsv; do
        prefixes="called_genes|calls|clinvar_calls|clinvar_lowmap_calls|clinvar_miscalls|clinvar_miscalls_lowmap|lowmap_called_genes|lowmap_calls|miscalled_clinvar_genes|miscalled_clinvar_genes_lowmap|miscalls|miscalls_lowmap|mixed_methyl_calls|pure_methyl_calls|resolved_calls|resolved_clinvar_genes|underfiltered_calls"
        metadata="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_([S0-9][^.]+)\\.([^.]+)\\.(bisulfite|emseq)\\.(grch38|t2t)\\.(bismark|bwameth)(\\.(filt|nofilt))?\\.tsv/\\4.\\5.\\6.\\7.\\8/")"
        processing_stage="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_([S0-9][^.]+).*.tsv/\\1/")"
        filt_status="$(echo "$tsv" | sed -E "s/counts_(($prefixes)(_mq0)?)_[S0-9][^.]+\\.[^.]+\\.[^.]+\\.[^.]+\\.[^.]+(\\.(filt|nofilt))?\\.tsv/\\5/")"
        output_metadata="$(echo "$tsv" | sed -E "s/counts_($prefixes)(_mq0)?_([S0-9][^.]+)\\.([^.]+)\\.(bisulfite|emseq)\\.(grch38|t2t)\\.(bismark|bwameth).*\\.tsv/\\3\t\\4\t\\5\t\\6\t\\7/")"
        if [ -z "$filt_status" ]; then 
            filt_status="combined";
        fi
        output_metadata="$filt_status\t$output_metadata"
        for context in CpG CHG CHH; do
            # appending a zero to the end of the list so no records returns 0, not empty
            total="$(grep "\tchr" "$tsv" | grep "\t${context}" | cat - <(echo "0") | cut -f 1 | paste -sd+ | bc)" 
            echo "$total\t${processing_stage}\t${context}\t${output_metadata}" >> "totals_${metadata}.tsv"
        done
    done <<< "!{counts_tsvs} "
    for tsv in totals_*.tsv; do
        ref_name="$(echo "$tsv" | cut -d "." -f 4)"
        c_count="$(grep "${ref_name}_all" !{c_counts} | cut -f 2)"        
        c_count_cpg="$(grep "${ref_name}_cpg" !{c_counts} | cut -f 2)"        
        awk -v "total=${c_count}" -v "total_cpg=${c_count_cpg}" 'BEGIN { OFS="\\t"; } { print $0, total, total_cpg }' "$tsv" > "$tsv".tmp;
        rm "$tsv";
        mv "$tsv".tmp "$tsv";
    done
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
    
    all_calls = calls.mix(clinvar_calls, clinvar_lowmap_calls, lowmap_calls, miscalls, pure_mixed_calls, underfiltered_calls)
    all_called_genes = called_genes.mix(lowmap_called_genes, miscalled_genes)

    refs = Channel.fromPath(params.refPath+"/*.fa")

    refs.view()

    count_per_chrom(all_calls)
    count_per_chrom_genes(all_called_genes)
    count_ref_cs(refs)
    combine_c_counts(count_ref_cs.out.collect()) 
    summarize_counts(count_per_chrom.out.mix(count_per_chrom_genes.out).collect(), combine_c_counts.out)
    sankey(summarize_counts.out.flatten())
}

