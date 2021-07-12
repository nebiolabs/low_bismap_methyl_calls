nextflow.enable.dsl=2

include { basicPipeline as GRCh38WGBSBismarkNoFilt } from './miscalls_one.nf'
include { basicPipeline as GRCh38WGBSBismarkFilt } from './miscalls_one.nf'

include { basicPipeline as GRCh38WGBSBwamethNoFilt } from './miscalls_one.nf'
include { basicPipeline as GRCh38WGBSBwamethFilt } from './miscalls_one.nf'

include { basicPipeline as GRCh38EMSeqBismarkNoFilt } from './miscalls_one.nf'
include { basicPipeline as GRCh38EMSeqBismarkFilt } from './miscalls_one.nf'

include { basicPipeline as GRCh38EMSeqBwamethNoFilt } from './miscalls_one.nf'
include { basicPipeline as GRCh38EMSeqBwamethFilt } from './miscalls_one.nf'

include { basicPipeline as T2TWGBSBismarkNoFilt } from './miscalls_one.nf'
include { basicPipeline as T2TWGBSBismarkFilt } from './miscalls_one.nf'

include { basicPipeline as T2TWGBSBwamethNoFilt } from './miscalls_one.nf'
include { basicPipeline as T2TWGBSBwamethFilt } from './miscalls_one.nf'

include { basicPipeline as T2TEMSeqBismarkNoFilt } from './miscalls_one.nf'
include { basicPipeline as T2TEMSeqBismarkFilt } from './miscalls_one.nf'

include { basicPipeline as T2TEMSeqBwamethNoFilt } from './miscalls_one.nf'
include { basicPipeline as T2TEMSeqBwamethFilt } from './miscalls_one.nf'

params.clinvarbedt2t="test_fixtures/chr18_clinvar_regions.bed"

params.bismapt2t="test_fixtures/chr18_bismap.bw"

params.bismapt2t_bbm="test_fixtures/chr18_bismap.bbm"

params.reft2t_fai="test_fixtures/chr18_test.fa.fai"

params.reft2t="test_fixtures/chr18_test.fa"

params.clinvarbed="test_fixtures/chr18_clinvar_regions.bed"

params.bismap="test_fixtures/chr18_bismap.bw"

params.bismap_bbm="test_fixtures/chr18_bismap.bbm"

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

params.tmpdir = "/tmp"

workflow {

  no_filter_mappability=params.no_filter_mappability

  bismap_bw_grch38 = file(params.bismap)
  bismap_bbm_grch38 = file(params.bismap_bbm)
  clinvar_regions_grch38 = file(params.clinvarbed)
  ref_fai_grch38 = file(params.ref_fai)
  ref_grch38 = file(params.ref)

  bismap_bw_t2t = file(params.bismapt2t)
  bismap_bbm_t2t = file(params.bismapt2t_bbm)
  clinvar_regions_t2t = file(params.clinvarbedt2t)
  ref_fai_t2t = file(params.reft2t_fai)
  ref_t2t = file(params.reft2t)

  GRCh38EMSeqBwamethNoFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, true, params.bams_bwameth_emseq, "grch38_bwameth_emseq_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  GRCh38EMSeqBwamethFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, false, params.bams_bwameth_emseq, "grch38_bwameth_emseq_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  GRCh38EMSeqBismarkNoFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, true, params.bams_bismark_emseq, "grch38_bismark_emseq_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  GRCh38EMSeqBismarkFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, false, params.bams_bismark_emseq, "grch38_bismark_emseq_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  GRCh38WGBSBwamethNoFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, true, params.bams_bwameth_wgbs, "grch38_bwameth_wgbs_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  GRCh38WGBSBwamethFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, false, params.bams_bwameth_wgbs, "grch38_bwameth_wgbs_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  GRCh38WGBSBismarkNoFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, true, params.bams_bismark_wgbs, "grch38_bismark_wgbs_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  GRCh38WGBSBismarkFilt(bismap_bbm_grch38, bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, false, params.bams_bismark_wgbs, "grch38_bismark_wgbs_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  T2TEMSeqBwamethNoFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bwameth_emseq, "t2t_bwameth_emseq_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  T2TEMSeqBwamethFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bwameth_emseq, "t2t_bwameth_emseq_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  T2TEMSeqBismarkNoFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bismark_emseq, "t2t_bismark_emseq_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  T2TEMSeqBismarkFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bismark_emseq, "t2t_bismark_emseq_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  T2TWGBSBwamethNoFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bwameth_wgbs, "t2t_bwameth_wgbs_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  T2TWGBSBwamethFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bwameth_wgbs, "t2t_bwameth_wgbs_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

  T2TWGBSBismarkNoFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, true, params.bamst2t_bismark_wgbs, "t2t_bismark_wgbs_nofilt", params.tmpdir, params.min_mapq, params.bismap_cutoff)
  T2TWGBSBismarkFilt(bismap_bbm_t2t, bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, false, params.bamst2t_bismark_wgbs, "t2t_bismark_wgbs_filt", params.tmpdir, params.min_mapq, params.bismap_cutoff)

}