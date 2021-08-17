nextflow.enable.dsl=2

params.bams="test_fixtures/foo.bam"
params.ref="test_fixtures/chr18_test.fa"
params.bismarkOnly=false
params.suffix="realigned"

process collateBams {

	conda 'samtools'

        penv 'smp'
        cpus 16

        input:
                each bam

        output:
                file "*_collated.bam"

        shell:
        '''
        samtools collate -@ 16 !{bam.replace("\\n", "")} !{bam.split("/")[-1].replace(".bam\\n", "").replace("\\n", "")}_collated
        '''
}

process extractReads {

	conda 'samtools'

        penv 'smp'
        cpus 16

        input:
                each bam

        output:
               tuple path("*_1.fastq"), path("*_2.fastq"), path("*_singletons.fastq"), path("*.bam")

        shell:
        '''
        samtools fastq -1 !{bam.getFileName().toString().split("/")[-1].replace(".bam", "")}_1.fastq -2 !{bam.getFileName().toString().split("/")[-1].replace(".bam", "")}_2.fastq -s !{bam.getFileName().toString().split("/")[-1].replace(".bam", "")}_singletons.fastq -@ 16 !{bam}; touch !{bam.getFileName().toString().split("/")[-1]};
        '''
}

process indexGenomeBwameth {

	conda 'bwameth samtools'


        input:
                file ref

        output:
               file 'genome_indexed.txt'
               tuple path("*.c2t"), path("*.c2t.amb"), path("*.c2t.ann"), path("*.c2t.bwt"), path("*.c2t.pac"), path("*.c2t.sa")

        shell:
        '''
        bwameth.py index !{ref}; echo "done" > genome_indexed.txt
        '''
}

process indexGenomeBismark {

	conda 'bismark bowtie2'


        input:
               path ref

        output:
               path 'Bisulfite_Genome'
               file 'genome_indexed.txt'

        shell:
        '''
        REALIGN_WORKDIR=`pwd`; echo $REALIGN_WORKDIR; bismark_genome_preparation $REALIGN_WORKDIR; echo "done" > genome_indexed.txt
        '''
}

process realignGenomeBismark {

        penv 'smp'
        cpus 16
        conda 'bismark bowtie2'


        publishDir 'bams', mode: 'copy'

        input:
                path ref
                tuple path('reads_1.fastq'), path('reads_2.fastq'), path('singletons.fastq'), path(bam)
                file 'genome_indexed.txt'
                path('Bisulfite_Genome')
                val suffix

        output:
               file "*_${suffix}.bam"

        shell:
        '''
        REALIGN_WORKDIR=`pwd`; echo $REALIGN_WORKDIR; bismark --parallel 3 --genome $REALIGN_WORKDIR -1 reads_1.fastq -2 reads_2.fastq; for f in `ls *.bam`; do mv $f !{bam.getFileName().toString().split("/")[-1].replace(".bam", "")}_bismark_!{suffix}.bam; done
        '''
}


process deduplicateBismark {

        penv 'smp'
        cpus 8
        conda 'bismark bowtie2'


        publishDir 'bams', mode: 'move'

        input:
                path aligned_bam
        output:
               file "*.deduplicated.bam"

        shell:
        '''
        deduplicate_bismark --bam !{aligned_bam}
        '''
}

process realignGenomeBwameth {

        penv 'smp'
        cpus 16
        conda 'bwameth samtools'


        publishDir 'bams', mode: 'move'

        input:
                path ref
                tuple path('reads_1.fastq'), path('reads_2.fastq'), path('singletons.fastq'), path(bam)
                file 'genome_indexed.txt'
                tuple path(c2t), path(c2t_amb), path(c2t_ann), path(c2t_bwt), path(c2t_pac), path(c2t_sa)
                val suffix

        output:
               file "*_${suffix}.bam"

        shell:
        '''
        bwameth.py --reference !{ref} -t 16 reads_1.fastq reads_2.fastq | samtools view -b - > !{bam.getFileName().toString().split("/")[-1].replace(".bam", "")}_bwameth_!{suffix}.bam
        '''
}

workflow {
  ref = file(params.ref)
  bam_names = params.bams.toString().split(",")

  println(bam_names)

  bam_paths = bam_names.collect({ "readlink -f ${it}".execute().text })
  bam_list = Channel.from(bam_paths)
  
  if( !params.bismarkOnly )
    indexGenomeBwameth(ref)
  indexGenomeBismark(ref)
  collateBams(bam_list)
  extractReads(collateBams.out)
  if( !params.bismarkOnly )
    realignGenomeBwameth(ref, extractReads.out, indexGenomeBwameth.out[0], indexGenomeBwameth.out[1], params.suffix)
  realignGenomeBismark(ref, extractReads.out, indexGenomeBismark.out[1], indexGenomeBismark.out[0], params.suffix)
  deduplicateBismark(realignGenomeBismark.out)
} 
