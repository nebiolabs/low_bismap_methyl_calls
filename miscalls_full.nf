nextflow.enable.dsl=2

include { basicPipeline } from './miscalls_one.nf'

params.clinvarbedt2t="test_fixtures/chr18_clinvar_regions.bed"

params.bismapt2t="test_fixtures/chr18_bismap.bw"

params.reft2t_fai="test_fixtures/chr18_test.fa.fai"

params.reft2t="test_fixtures/chr18_test.fa"

params.clinvarbed="test_fixtures/chr18_clinvar_regions.bed"

params.bismap="test_fixtures/chr18_bismap.bw"

params.ref_fai="test_fixtures/chr18_test.fa.fai"

params.ref="test_fixtures/chr18_test.fa"

params.bams_bwameth_wgbs="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamst2t_bwameth_wgbs="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bams_bismark_wgbs="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamst2t_bismark_wgbs="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bams_bwameth_emseq="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamst2t_bwameth_emseq="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bams_bismark_emseq="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamst2t_bismark_emseq="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.min_mapq=10

params.no_filter_mappability=false

params.bismap_cutoff=0.01


workflow {

  no_filter_mappability=params.no_filter_mappability

  bismap_bw_grch38 = file(params.bismap)
  clinvar_regions_grch38 = file(params.clinvarbed)
  ref_fai_grch38 = file(params.ref_fai)
  ref_grch38 = file(params.ref)

  bismap_bw_t2t = file(params.bismapt2t)
  clinvar_regions_t2t = file(params.clinvarbedt2t)
  ref_fai_t2t = file(params.reft2t_fai)
  ref_t2t = file(params.reft2t)

  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, true, params.bams_bwameth_emseq, "grch38_bwameth_emseq_nofilt")
  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, false, params.bams_bwameth_emseq, "grch38_bwameth_emseq_filt")

  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, true, params.bams_bismark_emseq, "grch38_bismark_emseq_nofilt")
  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, false, params.bams_bismark_emseq, "grch38_bismark_emseq_filt")

  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, true, params.bams_bwameth_wgbs, "grch38_bwameth_wgbs_nofilt")
  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, false, params.bams_bwameth_wgbs, "grch38_bwameth_wgbs_filt")

  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, true, params.bams_bismark_wgbs, "grch38_bismark_wgbs_nofilt")
  basicPipeline(bismap_bw, ref, ref_fai, clinvar_regions, false, params.bams_bismark_wgbs, "grch38_bismark_wgbs_filt")

  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bwameth_emseq, "t2t_bwameth_emseq_nofilt")
  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bwameth_emseq, "t2t_bwameth_emseq_filt")

  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bismark_emseq, "t2t_bismark_emseq_nofilt")
  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bismark_emseq, "t2t_bismark_emseq_filt")

  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bwameth_wgbs, "t2t_bwameth_wgbs_nofilt")
  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bwameth_wgbs, "t2t_bwameth_wgbs_filt")

  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bismark_wgbs, "t2t_bismark_wgbs_nofilt")
  basicPipeline(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bismark_wgbs, "t2t_bismark_wgbs_filt")

}