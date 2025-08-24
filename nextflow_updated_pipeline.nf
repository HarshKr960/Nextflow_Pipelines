nextflow.enable.dsl=2

params.reads     = '/home/Downloads/HARSH_TEST_AREA/git_test/data/second_data/*_{1,2}.fastq.gz'
params.adapters  = '/home/miniconda3/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa'
params.index     = '/home/ref/grch38/genome'  // HISAT2 index prefix
params.gtf       = '/home/ref/Homo_sapiens.GRCh38.113.gtf'
params.output    = '/home/Downloads/HARSH_TEST_AREA/git_test/data/second_data/nextflow_output'
params.threads   = 16

process SetupDirectories {
    output:
    path "setup_complete.txt"

    script:
    """
    mkdir -p ${params.output}/{fastqc,trimmed,hisat2,bam,counts}
    echo "Setup completed at \$(date)" > setup_complete.txt
    """
}

process FastQC {
    publishDir "${params.output}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.{html,zip}")

    script:
    """
    fastqc --threads ${params.threads} -o . ${reads}
    """
}

process Trimmomatic {
    publishDir "${params.output}/trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_paired.fq.gz")

    script:
    """
    read1=${reads[0]}
    read2=${reads[1]}
    base=${sample_id}

    trimmomatic PE -threads ${params.threads} \\
        \$read1 \$read2 \\
        \${base}_1_paired.fq.gz \${base}_1_unpaired.fq.gz \\
        \${base}_2_paired.fq.gz \${base}_2_unpaired.fq.gz \\
        ILLUMINACLIP:${params.adapters}:2:30:10 \\
        SLIDINGWINDOW:4:20 \\
        MINLEN:36
    """
}

process HISAT2 {
    publishDir "${params.output}/hisat2", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}.hisat2_mapstats.txt")


    script:
    """
    read1=${trimmed_reads[0]}
    read2=${trimmed_reads[1]}

    hisat2 -p ${params.threads} \\
        --very-sensitive \\
        -x ${params.index} \\
        -1 \$read1 -2 \$read2 \\
        --rna-strandness RF \\
        --pen-noncansplice 14 \\
        --mp 1,0 \\
        --sp 3,1 \\
        --new-summary \\
        -S ${sample_id}.sam \\
        2> ${sample_id}.hisat2_mapstats.txt
    """
}

process SamtoolsSort {
    publishDir "${params.output}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file), path(mapstats)


    output:
    tuple val(sample_id), path("*.sorted.bam")

    script:
    """
    samtools view -@ ${params.threads} -bS $sam_file | \\
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam
    """
}

process FeatureCounts {
    publishDir "${params.output}/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "featureCounts.txt"
    path "featureCounts.txt.summary"

    script:
    """
    featureCounts -T ${params.threads} \\
        -a ${params.gtf} \\
        -o featureCounts.txt \\
        -p --countReadPairs -B -C \\
        -t exon -g gene_id \\
        ${bam_files.join(' ')}
    """
}

workflow {
    setup = SetupDirectories()

    reads_ch = Channel.fromFilePairs(params.reads, size: 2)

    fastqc_results = FastQC(reads_ch)

    trimmed_reads = Trimmomatic(reads_ch)

    aligned_reads = HISAT2(trimmed_reads)

    sorted_bams = SamtoolsSort(aligned_reads)

    bam_ch = sorted_bams.map { sample_id, bam -> bam }

    FeatureCounts(bam_ch.collect())
}
