nextflow.enable.dsl=2

params.bismark = "true"
params.bwameth = "true"

process indexGenomeBwameth {

        penv 'smp'
        cpus 16
        conda 'bwameth==0.2.7 samtools==1.19'

        storeDir "${params.data_path}/index_bwameth/${suffix}"

        input:
                tuple path(ref), val(suffix)

        output:
		tuple path(ref), path("*.c2t"), path("*.c2t.amb"), path("*.c2t.ann"), path("*.c2t.bwt"), path("*.c2t.pac"), path("*.c2t.sa"), val(suffix)

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
		bwameth.py index !{ref};
		'''
}

process indexGenomeBismark {

        penv 'smp'
        cpus 16
	conda 'bismark==0.24.1 bowtie2==2.5.2'

        storeDir "${params.data_path}/index_bismark/${suffix}"

        input:
                tuple path(ref), val(suffix)

        output:
		tuple path(ref), path('Bisulfite*'), val(suffix)
        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
		bismark_genome_preparation "$(pwd)" --parallel 8;
                '''
}

process realignGenomeBismark {

        penv 'smp'
        cpus 16
        conda 'bismark==0.24.1 bowtie2==2.5.2 samtools==1.19 gzip==1.13 seqtk==1.4 fastp==0.20.1'

        input:
                tuple path(reads_1), path(reads_2), path(ref), path('Bisulfite_Genome'), val(suffix)

        output:
		tuple file("*_${suffix}.bam"), val(suffix)

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
                cleanup() {
                    rm reads_1.fastq reads_2.fastq 
                }
                trap cleanup EXIT
                # reads_1 is actually a BAM and reads_2 is actually a BAI
                if [[ "!{reads_1}" = *.bam ]]; then
                    reads1_cmd="samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -2 /dev/null !{reads_1}"
                    reads2_cmd="samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -1 /dev/null !{reads_1}"
                else
                    reads1_cmd="zcat -f \"!{reads_1}\""
                    reads2_cmd="zcat -f \"!{reads_2}\""
                fi
                fastp -l 2 -Q --in1 <($reads1_cmd) --in2 <($reads2_cmd) --out1 reads_1.fastq --out2 reads_2.fastq --overrepresentation_analysis 2> fastp.stderr
		bismark --parallel 3 --temp_dir "$TMPDIR" "." -1 reads_1.fastq -2 reads_2.fastq; for f in `ls *_bismark*.bam`; do mv $f !{reads_1.getFileName().toString().split("/")[-1].replace(".bam", "").replace(".fastq", "").replace(".fq", "").replace(".gz", "")}_bismark_!{suffix}.bam; done
                '''
}

process deduplicateBismark {

        penv 'smp'
        cpus 8
        conda 'bismark==0.24.1 bowtie2==2.5.2'

        input:
                tuple path(aligned_bam), val(suffix)
        output:
		tuple file("*.deduplicated.bam"), val(suffix)

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
                bam_name_noext="$(basename -s ".bam" !{aligned_bam})"
		deduplicate_bismark -p --bam !{aligned_bam}
		'''
}

process realignGenomeBwameth {

        penv 'smp'
        cpus 16
        conda 'bwameth==0.2.7 samtools==1.19 gzip==1.13 seqtk==1.4 fastp==0.20.1'

        input:
		tuple path(reads_1), path(reads_2), path(ref), path(c2t), path(c2t_amb), path(c2t_ann), path(c2t_bwt), path(c2t_pac), path(c2t_sa), val(suffix)

        output:
		tuple file("*_${suffix}.bam"), val(suffix)

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
                cleanup() {
                    rm reads_1.fastq reads_2.fastq 
                }
                trap cleanup EXIT
                #mkfifo reads1_trimmed reads2_trimmed
                # reads_1 is actually a BAM and reads_2 is actually a BAI
                if [[ "!{reads_1}" = *.bam ]]; then
                    reads1_cmd="samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -2 /dev/null !{reads_1}"
                    reads2_cmd="samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -1 /dev/null !{reads_1}"
                else
                    reads1_cmd="zcat -f \"!{reads_1}\""
                    reads2_cmd="zcat -f \"!{reads_2}\""
                fi
                fastp -l 2 -Q --in1 <($reads1_cmd) --in2 <($reads2_cmd) --out1 reads_1.fastq --out2 reads_2.fastq --overrepresentation_analysis 2> fastp.stderr
		bwameth.py -p --reference !{ref} -t 16 reads_1.fastq reads_2.fastq | samtools view -b - > !{reads_1.getFileName().toString().split("/")[-1].replace(".bam", "").replace(".fastq", "").replace(".fq", "").replace(".gz", "")}_bwameth_!{suffix}.bam
		'''
}

process addReadGroup {

        penv 'smp'
        cpus 1
        conda "samtools==1.19"

        input:
                tuple file(input_bam), val(suffix)

        output:
                tuple file("*.rg.bam"), val(suffix)

        shell:
                '''
                mkdir tmp
                export TMPDIR="$(pwd)/tmp"
                samtools addreplacerg -r "@RG\tID:-\tSM:-" -o "$(basename -s .bam "!{input_bam}").rg.bam" "!{input_bam}"
                '''
}

process indexBismark {

        penv 'smp'
        cpus 8
        conda 'samtools==1.19'

        publishDir "${params.output_path}/bams_${suffix}_bismark", mode: 'copy'

        input:
                tuple path(aligned_bam), val(suffix)
        output:
                tuple file("*.sorted.bam"), file("*.sorted.bam.bai")

        shell:
                '''
                mkdir tmp
                export TMPDIR="$(pwd)/tmp"
                bam_name_noext="$(basename -s ".bam" !{aligned_bam})"
                samtools sort -@ 8 !{aligned_bam} > "${bam_name_noext}.sorted.bam"
                samtools index -@ 8 "${bam_name_noext}.sorted.bam"
                '''
}

process sortBwameth {

        penv 'smp'
        cpus 8
        conda 'samtools==1.19'

        input:
                tuple file(aligned_bam), val(suffix)
        output:
                tuple file("*.sorted.bam"), val(suffix)

        shell:
                '''
                mkdir tmp
                export TMPDIR="$(pwd)/tmp"
                bam_name_noext="$(basename -s ".bam" !{aligned_bam})"
                samtools sort -@ 8 !{aligned_bam} > "${bam_name_noext}.sorted.bam"
                '''
}


process markDupsBwameth {

        penv 'smp'
        cpus 8
        conda 'picard==2.18.29 samtools==1.19'

        publishDir "${params.output_path}/bams_${suffix}_bwameth", mode: 'copy'

        input:
                tuple file(sorted_bam), val(suffix)
        output:
                tuple file("*.md.sorted.bam"), file("*.md.sorted.bam.bai")

        shell:
                '''
                mkdir tmp
                export TMPDIR="$(pwd)/tmp"
                bam_name_noext="$(basename -s ".sorted.bam" !{sorted_bam})"
                picard MarkDuplicates \
                I=!{sorted_bam} \
                O="${bam_name_noext}.md.sorted.bam" \
                M="${bam_name_noext}_dup_metrics.txt";
                samtools index -@ 8 "${bam_name_noext}.md.sorted.bam"
                '''
}

workflow {
 
  fastqs = Channel.fromFilePairs("${params.input_path}/*/*_{1,2}.f*q.gz").map { it[1] }
  bams = Channel.fromPath("${params.input_path}/*/*.bam").map { [it, it+".bai"] }
  seqs = fastqs.mix(bams)

  refs_list = params.refs.split(",")
  suffixes_list = params.suffixes.split(",")
 
  refs_count = refs_list.size()
  ref_pairs = []

  for(int i=0; i<refs_count; i++)
  {
      ref_pairs.add([file(refs_list[i]).toAbsolutePath(), suffixes_list[i]])
  }

  ref_pairs_ch = Channel.fromList(ref_pairs)

  if(params.bwameth == "true")
  {
      //index the genome
      indexGenomeBwameth(ref_pairs_ch)
      //combine channels
      fastqs_bwameth = seqs.combine(indexGenomeBwameth.out)
      //realign
      realignGenomeBwameth(fastqs_bwameth)
      //sort
      sortBwameth(realignGenomeBwameth.out)
      //deduplicate and index
      markDupsBwameth(sortBwameth.out)
  }
  if(params.bismark == "true")
  {
      //index the genome
      indexGenomeBismark(ref_pairs_ch)
      //combine channels
      fastqs_bismark = seqs.combine(indexGenomeBismark.out)
      //realign
      realignGenomeBismark(fastqs_bismark)
      //add @RG tag
      addReadGroup(realignGenomeBismark.out)
      //deduplicate
      deduplicateBismark(addReadGroup.out)
      //sort and index
      indexBismark(deduplicateBismark.out.mix(addReadGroup.out))
  } 
} 
