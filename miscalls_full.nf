nextflow.enable.dsl=2

include { refPreproc as PreprocRefGRCh38 } from './miscalls_one.nf'
include { refPreproc as PreprocRefT2T } from './miscalls_one.nf'

include { runPair as GRCh38WGBSBismark } from './miscalls_pair.nf'
include { runPair as GRCh38WGBSBwameth } from './miscalls_pair.nf'
include { runPair as GRCh38EMSeqBismark } from './miscalls_pair.nf'
include { runPair as GRCh38EMSeqBwameth } from './miscalls_pair.nf'
include { runPair as T2TWGBSBismark } from './miscalls_pair.nf'
include { runPair as T2TWGBSBwameth } from './miscalls_pair.nf'
include { runPair as T2TEMSeqBismark} from './miscalls_pair.nf'
include { runPair as T2TEMSeqBwameth } from './miscalls_pair.nf'

params.clinvarBedT2t="test_fixtures/chr18_clinvar_regions.bed"

params.bismapBwT2t="test_fixtures/chr18_bismap.bw"

params.bismapBbmT2t="test_fixtures/chr18_bismap.bbm"

params.refGrch38FaiT2t="test_fixtures/chr18_test.fa.fai"

params.refT2t="test_fixtures/chr18_test.fa"

params.clinvarBedGrch38="test_fixtures/chr18_clinvar_regions.bed"

params.bismapBwGrch38="test_fixtures/chr18_bismap.bw"

params.bismapBbmGrch38="test_fixtures/chr18_bismap.bbm"

params.refGrch38FaiGrch38="test_fixtures/chr18_test.fa.fai"

params.refGrch38="test_fixtures/chr18_test.fa"

params.bamsBwamethWgbsGrch38="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBwamethWgbsT2t="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBismarkWgbsGrch38="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBismarkWgbsT2t="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBwamethEmseqGrch38="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBwamethEmseqT2t="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBismarkEmseqGrch38="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.bamsBismarkEmseqT2t="test_fixtures/chr18_bwameth_fixed.bam,test_fixtures/chr18_bismark_fixed.bam"

params.minMapq=10

params.noFilterMappability=false

params.bismapCutoff=0.01

params.tmpdir = "/tmp"

params.massMapping = "test_fixtures/mass_mapping.csv"

params.copy = false

workflow {

    no_filter_mappability=params.noFilterMappability

    bismap_bw_grch38 = file(params.bismapBwGrch38)
    bismap_bbm_grch38 = file(params.bismapBbmGrch38)
    clinvar_regions_grch38 = file(params.clinvarBedGrch38)
    ref_fai_grch38 = file(params.refGrch38FaiGrch38)
    ref_grch38 = file(params.refGrch38)

    bismap_bw_t2t = file(params.bismapBwT2t)
    bismap_bbm_t2t = file(params.bismapBbmT2t)
    clinvar_regions_t2t = file(params.clinvarBedT2t)
    ref_fai_t2t = file(params.refGrch38FaiT2t)
    ref_t2t = file(params.refT2t)

    mass_mapping = file(params.massMapping)

    publishMode = "symlink"

    if(params.copy)
    {
        publishMode = "copy"
    }

    PreprocRefGRCh38(bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38, params.tmpdir, params.minMapq, params.bismapCutoff)

    PreprocRefT2T(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, params.tmpdir, params.minMapq, params.bismapCutoff)

    GRCh38EMSeqBwameth(bismap_bbm_grch38, ref_grch38, ref_fai_grch38, params.bamsBwamethEmseqGrch38, "grch38", "bwameth", "emseq", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefGRCh38.out[0], PreprocRefGRCh38.out[1], PreprocRefGRCh38.out[2], publishMode)

    GRCh38EMSeqBismark(bismap_bbm_grch38, ref_grch38, ref_fai_grch38, params.bamsBismarkEmseqGrch38, "grch38", "bismark", "emseq", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefGRCh38.out[0], PreprocRefGRCh38.out[1], PreprocRefGRCh38.out[2], publishMode)

    GRCh38WGBSBwameth(bismap_bbm_grch38, ref_grch38, ref_fai_grch38, params.bamsBwamethWgbsGrch38, "grch38", "bwameth", "wgbs", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefGRCh38.out[0], PreprocRefGRCh38.out[1], PreprocRefGRCh38.out[2], publishMode)

    GRCh38WGBSBismark(bismap_bbm_grch38, ref_grch38, ref_fai_grch38, params.bamsBismarkWgbsGrch38, "grch38", "bismark", "wgbs", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefGRCh38.out[0], PreprocRefGRCh38.out[1], PreprocRefGRCh38.out[2], publishMode)

    T2TEMSeqBwameth(bismap_bbm_t2t, ref_t2t, ref_fai_t2t, params.bamsBwamethEmseqT2t, "t2t", "bwameth", "emseq", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefT2T.out[0], PreprocRefT2T.out[1], PreprocRefT2T.out[2], publishMode)

    T2TEMSeqBismark(bismap_bbm_t2t, ref_t2t, ref_fai_t2t, params.bamsBismarkEmseqT2t, "t2t", "bismark", "emseq", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefT2T.out[0], PreprocRefT2T.out[1], PreprocRefT2T.out[2], publishMode)

    T2TWGBSBwameth(bismap_bbm_t2t, ref_t2t, ref_fai_t2t, params.bamsBwamethWgbsT2t, "t2t", "bwameth", "wgbs", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefT2T.out[0], PreprocRefT2T.out[1], PreprocRefT2T.out[2], publishMode)

    T2TWGBSBismark(bismap_bbm_t2t, ref_t2t, ref_fai_t2t, params.bamsBismarkWgbsT2t, "t2t", "bismark", "wgbs", params.tmpdir, params.minMapq, params.bismapCutoff, mass_mapping, PreprocRefT2T.out[0], PreprocRefT2T.out[1], PreprocRefT2T.out[2], publishMode)

}