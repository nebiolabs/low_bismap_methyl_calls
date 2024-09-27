nextflow.enable.dsl=2

process indexGenomeBwameth {

        penv 'smp'
        cpus 16
        conda 'bwameth=0.2.7 samtools=1.18'

        storeDir "${params.data_path}/index_bwameth"

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
	conda 'bismark=0.24.1 bowtie2=2.5.2'

        storeDir "${params.data_path}/index_bismark"

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
        conda 'bismark=0.24.1 bowtie2=2.5.2 gzip=1.13 seqtk=1.4'


        publishDir "${params.data_path}/bams_${suffix}_bismark", mode: 'copy'

        input:
                tuple path(reads_1), path(reads_2), path(ref), path('Bisulfite_Genome'), val(suffix)

        output:
		tuple file("*_${suffix}.bam"), val(suffix)

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
                # reads_1 is actually a BAM and reads_2 is actually a BAI
                if [[ "!{reads_1}" = *.bam ]]; then
                    samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -1 reads_1.fastq -2 reads_2.fastq !{reads_1} 
                else
                    seqtk mergepe <(zcat -f "!{reads_1}") <(zcat -f "!{reads_2}") | tee >(seqtk seq -1 > reads_1.fastq) | seqtk seq -2 > reads_2.fastq
                fi
		bismark --parallel 3 --temp_dir "$TMPDIR" "." -1 reads_1.fastq -2 reads_2.fastq; for f in `ls *_bismark*.bam`; do mv $f !{reads_1.getFileName().toString().split("/")[-1].replace(".bam", "").replace(".fastq", "").replace(".fq", "")}_bismark_!{suffix}.bam; done
                rm reads_1.fastq
                rm reads_2.fastq
                '''
}

process deduplicateBismark {

        penv 'smp'
        cpus 8
        conda 'bismark=0.24.1 bowtie2=2.5.2'


        publishDir "${params.data_path}/bams_${suffix}_bismark", mode: 'copy'

        input:
                tuple path(aligned_bam), val(suffix)
        output:
		file "*.deduplicated.bam"

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
        conda 'bwameth=0.2.7 samtools=1.18 gzip=1.13 seqtk=1.4'


        publishDir "${params.data_path}/bams_${suffix}_bwameth", mode: 'copy'

        input:
		tuple path(reads_1), path(reads_2), path(ref), path(c2t), path(c2t_amb), path(c2t_ann), path(c2t_bwt), path(c2t_pac), path(c2t_sa), val(suffix)

        output:
		file "*_${suffix}.bam"

        shell:
		'''
                mkdir tmp
		export TMPDIR="$(pwd)/tmp"
                # reads_1 is actually a BAM and reads_2 is actually a BAI
                if [[ "!{reads_1}" = *.bam ]]; then
                    samtools collate -@ 2 -u -O !{reads_1} | samtools fastq -s /dev/null -0 /dev/null -1 reads_1.fastq -2 reads_2.fastq !{reads_1}
		    bwameth.py -p --reference !{ref} -t 16 reads_1.fastq reads_2.fastq | samtools view -b - > !{reads_1.getFileName().toString().split("/")[-1].replace(".bam", "").replace(".fastq", "").replace(".fq", "")}_bwameth_!{suffix}.bam
                    rm reads_1.fastq
                    rm reads_2.fastq
                else
		    bwameth.py --reference !{ref} -t 16 <(zcat -f !{reads_1}) <(zcat -f !{reads_2}) | samtools view -b - > !{reads_1.getFileName().toString().split("/")[-1].replace(".bam", "").replace(".fastq", "").replace(".fq", "")}_bwameth_!{suffix}.bam
                fi
		'''
}

process indexBismark {

        penv 'smp'
        cpus 8
        conda 'samtools==1.6'

        publishDir "${params.data_path}/bams_${suffix}_bismark", mode: 'copy'

        input:
                path aligned_bam
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
        conda 'samtools==1.6'

        input:
                file  aligned_bam
        output:
                file "*.sorted.bam"

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

        publishDir "${params.data_path}/bams_${suffix}_bwameth", mode: 'copy'

        input:
                file sorted_bam
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

  //index the genome
  //indexGenomeBwameth(ref_pairs_ch)
  indexGenomeBismark(ref_pairs_ch)
 
  //combine channels
  //fastqs_bwameth = seqs.combine(indexGenomeBwameth.out)
  fastqs_bismark = seqs.combine(indexGenomeBismark.out)

  //realign
  //realignGenomeBwameth(fastqs_bwameth)
  realignGenomeBismark(fastqs_bismark)
  deduplicateBismark(realignGenomeBismark.out)
  indexBismark(deduplicateBismark.out)
  //sortBwameth(realignGenomeBwameth.out)
  //markDupsBwameth(sortBwameth.out)

} 
