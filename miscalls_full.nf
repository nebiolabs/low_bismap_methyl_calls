nextflow.enable.dsl=2

include { refPreproc as PreprocRefGRCh38 } from './miscalls_one.nf'
include { refPreproc as PreprocRefT2T } from './miscalls_one.nf'

include { runPair as GRCh38 } from './miscalls_pair.nf'
include { runPair as T2T } from './miscalls_pair.nf'

params.clinvarBedT2t="test_fixtures/chr18_clinvar_regions.bed"

params.bismapBwT2t="test_fixtures/chr18_bismap.bw"

params.bismapBbmT2t="test_fixtures/chr18_bismap.bbm"

params.refFaiT2t="test_fixtures/chr18_test.fa.fai"

params.refT2t="test_fixtures/chr18_test.fa"

params.clinvarBedGrch38="test_fixtures/chr18_clinvar_regions.bed"

params.bismapBwGrch38="test_fixtures/chr18_bismap.bw"

params.bismapBbmGrch38="test_fixtures/chr18_bismap.bbm"

params.refFaiGrch38="test_fixtures/chr18_test.fa.fai"

params.refGrch38="test_fixtures/chr18_test.fa"

params.bamsDir="/this/does/not/exist"

params.minMapq=10

params.noFilterMappability=false

params.bismapCutoff=0.01

params.tmpdir = "/tmp"

params.hardlink = false

params.outputDir="methyl_calls"

workflow {

    no_filter_mappability=params.noFilterMappability

    bismap_bw_grch38 = file(params.bismapBwGrch38)
    bismap_bbm_grch38 = file(params.bismapBbmGrch38)
    clinvar_regions_grch38 = file(params.clinvarBedGrch38)
    ref_fai_grch38 = file(params.refFaiGrch38)
    ref_grch38 = file(params.refGrch38)

    bismap_bw_t2t = file(params.bismapBwT2t)
    bismap_bbm_t2t = file(params.bismapBbmT2t)
    clinvar_regions_t2t = file(params.clinvarBedT2t)
    ref_fai_t2t = file(params.refFaiT2t)
    ref_t2t = file(params.refT2t)

    publishMode = "symlink"

    if(params.hardlink)
    {
        publishMode = "link"
    }

    bam_pairs_grch38 = Channel.fromFilePairs(params.bamsDir+"*.*.*.grch38.*.{bam,bam.bai}").map { seqName, files -> [seqName, files[0], files[1]] }
    bam_pairs_t2t = Channel.fromFilePairs(params.bamsDir+"*.*.*.t2t.*.{bam,bam.bai}").map { seqName, files -> [seqName, files[0], files[1]] }


    bam_pairs_grch38.view()

    PreprocRefGRCh38(bismap_bw_grch38, ref_grch38, ref_fai_grch38, clinvar_regions_grch38,
                     params.tmpdir, params.minMapq, params.bismapCutoff)
    grch38_clinvar_bismap_low = PreprocRefGRCh38.out[0]
    grch38_combined_low_regions = PreprocRefGRCh38.out[1]
    grch38_stats = PreprocRefGRCh38.out[2]

    PreprocRefT2T(bismap_bw_t2t, ref_t2t, ref_fai_t2t, clinvar_regions_t2t, 
                  params.tmpdir, params.minMapq, params.bismapCutoff)
    t2t_clinvar_bismap_low = PreprocRefT2T.out[0]
    t2t_combined_low_regions = PreprocRefT2T.out[1]
    t2t_stats = PreprocRefT2T.out[2]

    GRCh38(bismap_bbm_grch38, ref_grch38, ref_fai_grch38,
                       bam_pairs_grch38,
                       params.tmpdir, params.minMapq, params.bismapCutoff,
                       grch38_clinvar_bismap_low,
                       grch38_combined_low_regions, grch38_stats, publishMode)

    T2T(bismap_bbm_t2t, ref_t2t, ref_fai_t2t,
                    bam_pairs_t2t,
                    params.tmpdir, params.minMapq, params.bismapCutoff,
                    t2t_clinvar_bismap_low,
                    t2t_combined_low_regions, t2t_stats, publishMode)

}
