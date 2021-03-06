configfile: "config.yaml"

import os
from os.path import basename

OUTPUT_DIR = config['output_dir']
FASTQ_DIR = os.path.join(OUTPUT_DIR, config['fastq_dir'])
GENOME_DIR = os.path.join(OUTPUT_DIR, config['genome_dir'])
CLASS_DIR = os.path.join(OUTPUT_DIR, config['class_dir'])
INSERTSEQ_DIR = os.path.join(OUTPUT_DIR, config['insertseq_dir'])

GENOME_ALIGNMENT_DIR = os.path.join(OUTPUT_DIR, config['genome_alignment_dir'])
INSERTSEQ_ALIGNMENT_DIR = os.path.join(OUTPUT_DIR, config['insertseq_alignment_dir'])
GENOME_INSERTSIZE_STATS = os.path.join(OUTPUT_DIR, config['genome_insertsize_stats_dir'])

GENOME_FLANKS_DIR = os.path.join(OUTPUT_DIR, config['genome_flanks_dir'])
INSERTSEQ_FLANKS_DIR = os.path.join(OUTPUT_DIR, config['insertseq_flanks_dir'])

WC_fastqs = glob_wildcards(os.path.join(FASTQ_DIR, "{sample}.{pair}.fq.gz"))
WC_genomes = glob_wildcards(os.path.join(GENOME_DIR, "{genome}.fasta"))
WC_inseqs = glob_wildcards(os.path.join(INSERTSEQ_DIR, "{insertseq}.fasta"))

SAMPLES = set(WC_fastqs.sample)

GENOMES = WC_genomes.genome
INSERTSEQS = WC_inseqs.insertseq

print(SAMPLES)
print(GENOMES)
print(INSERTSEQS)

rule all:
    input:
        expand("%s/{sample}.{genome}.insert_size_metrics.tsv" % GENOME_INSERTSIZE_STATS, sample=SAMPLES, genome=GENOMES),
        expand("%s/{sample}.{genome}.{insertseq}.flanks.bam.bai" % GENOME_FLANKS_DIR, sample=SAMPLES, genome=GENOMES, insertseq=INSERTSEQS),
        expand("%s/{sample}.{genome}.{insertseq}.noflanks.bam.bai" % GENOME_FLANKS_DIR, sample=SAMPLES, genome=GENOMES, insertseq=INSERTSEQS),
        expand("%s/{sample}.{genome}.{insertseq}.flanks.bam.bai" % INSERTSEQ_FLANKS_DIR, sample=SAMPLES, genome=GENOMES, insertseq=INSERTSEQS),
        expand("%s/{sample}.{genome}.{insertseq}.noflanks.bam.bai" % INSERTSEQ_FLANKS_DIR, sample=SAMPLES, genome=GENOMES, insertseq=INSERTSEQS),
        #expand("%s/{sample}.{genome}.{insertseq}.insert_size_metrics.txt" % GENOME_INSERTSIZE_STATS, sample=SAMPLES, genome=GENOMES, insertseq=INSERTSEQS)
        #expand("%s/{sample}.{genome}.complete" % GENOME_ALIGNMENT_DIR, sample=SAMPLES, genome=GENOMES),

rule bowtie_build_genomes:
    input:
        "%s/{genome}.fasta" % GENOME_DIR
    output:
        "%s/{genome}.fasta.1.bt2" % GENOME_DIR,
        "%s/{genome}.fasta.2.bt2" % GENOME_DIR,
        "%s/{genome}.fasta.3.bt2" % GENOME_DIR,
        "%s/{genome}.fasta.4.bt2" % GENOME_DIR,
        "%s/{genome}.fasta.rev.1.bt2" % GENOME_DIR,
        "%s/{genome}.fasta.rev.2.bt2" % GENOME_DIR
    shell:
        "bowtie2-build {input} {input}"
        
        
rule bowtie_build_insertseqs:
    input:
        "%s/{insertseq}.fasta" % INSERTSEQ_DIR
    output:
        "%s/{insertseq}.fasta.1.bt2" % INSERTSEQ_DIR,
        "%s/{insertseq}.fasta.2.bt2" % INSERTSEQ_DIR,
        "%s/{insertseq}.fasta.3.bt2" % INSERTSEQ_DIR,
        "%s/{insertseq}.fasta.4.bt2" % INSERTSEQ_DIR,
        "%s/{insertseq}.fasta.rev.1.bt2" % INSERTSEQ_DIR,
        "%s/{insertseq}.fasta.rev.2.bt2" % INSERTSEQ_DIR
    shell:
        "bowtie2-build {input} {input}"

rule bowtie_align_genome:
    input:
        fq1="%s/{sample}.R1.fq.gz" % FASTQ_DIR,
        fq2="%s/{sample}.R2.fq.gz" % FASTQ_DIR,
        genome="%s/{genome}.fasta" % GENOME_DIR,
        in1="%s/{genome}.fasta.1.bt2" % GENOME_DIR,
        in2="%s/{genome}.fasta.2.bt2" % GENOME_DIR,
        in3="%s/{genome}.fasta.3.bt2" % GENOME_DIR,
        in4="%s/{genome}.fasta.4.bt2" % GENOME_DIR,
        in5="%s/{genome}.fasta.rev.1.bt2" % GENOME_DIR,
        in6="%s/{genome}.fasta.rev.2.bt2" % GENOME_DIR
    output:
        sam = "%s/{sample}.{genome}.sam" % GENOME_ALIGNMENT_DIR,
        met = "%s/{sample}.{genome}.metrics.txt" % GENOME_ALIGNMENT_DIR,
        complete = "%s/{sample}.{genome}.complete" % GENOME_ALIGNMENT_DIR
    threads:
        config['bowtie2_threads']
    shell:
        "bowtie2 -x {input.genome} -1 {input.fq1} -2 {input.fq2} --met-file {output.met} -X 5000 "
        "-S {output.sam} --local --quiet --reorder -p {threads}; "
        "grep -v -P '^[^\t]+\t[^\t]+\t\*' {output.sam} >> {output.sam}.tmp; "
        "mv {output.sam}.tmp {output.sam}; "
        "echo 'COMPLETE' > {output.complete}"

rule bowtie_align_insertseq:
    input:
        fq1="%s/{sample}.R1.fq.gz" % FASTQ_DIR,
        fq2="%s/{sample}.R2.fq.gz" % FASTQ_DIR,
        genome="%s/{insertseq}.fasta" % INSERTSEQ_DIR,
        in1="%s/{insertseq}.fasta.1.bt2" % INSERTSEQ_DIR,
        in2="%s/{insertseq}.fasta.2.bt2" % INSERTSEQ_DIR,
        in3="%s/{insertseq}.fasta.3.bt2" % INSERTSEQ_DIR,
        in4="%s/{insertseq}.fasta.4.bt2" % INSERTSEQ_DIR,
        in5="%s/{insertseq}.fasta.rev.1.bt2" % INSERTSEQ_DIR,
        in6="%s/{insertseq}.fasta.rev.2.bt2" % INSERTSEQ_DIR
    output:
        sam = "%s/{sample}.{insertseq}.sam" % INSERTSEQ_ALIGNMENT_DIR,
        met = "%s/{sample}.{insertseq}.metrics.txt" % INSERTSEQ_ALIGNMENT_DIR,
        complete = "%s/{sample}.{insertseq}.complete" % INSERTSEQ_ALIGNMENT_DIR
    threads:
        config['bowtie2_threads']
    shell:
        "bowtie2 -x {input.genome} -1 {input.fq1} -2 {input.fq2} --met-file {output.met} -X 5000 "
        "-S {output.sam} --local --quiet --reorder -p {threads}; "
        "grep -v -P '^[^\t]+\t[^\t]+\t\*' {output.sam} > {output.sam}.tmp; "
        "mv {output.sam}.tmp {output.sam}; "
        "echo 'COMPLETE' > {output.complete}"

rule get_genome_insert_size_stats:
    input:
        sam = "%s/{sample}.{genome}.sam" % GENOME_ALIGNMENT_DIR
    output:
        stats = "%s/{sample}.{genome}.insert_size_metrics.tsv" % GENOME_INSERTSIZE_STATS,
        counts = "%s/{sample}.{genome}.insert_size_counts.tsv" % GENOME_INSERTSIZE_STATS
    params:
        sample='{sample}'
    shell:
        "python scripts/insert_size_stats.py {input.sam} {output.stats} {output.counts} {params.sample}"

rule get_insertseq_flanks:
    input:
        class1 = "%s/{sample}.R1.class.gz" % CLASS_DIR,
        class2 = "%s/{sample}.R2.class.gz" % CLASS_DIR,
        genome_sam = "%s/{sample}.{genome}.sam" % GENOME_ALIGNMENT_DIR,
        genome_met = "%s/{sample}.{genome}.metrics.txt" % GENOME_ALIGNMENT_DIR,
        insertseq_sam = "%s/{sample}.{insertseq}.sam" % INSERTSEQ_ALIGNMENT_DIR,
        insertseq_met = "%s/{sample}.{insertseq}.metrics.txt" % INSERTSEQ_ALIGNMENT_DIR,
        complete_genome = "%s/{sample}.{genome}.complete" % GENOME_ALIGNMENT_DIR,
        complete_insertseq = "%s/{sample}.{insertseq}.complete" % INSERTSEQ_ALIGNMENT_DIR
    output:
        genome_flanks_bam = "%s/{sample}.{genome}.{insertseq}.flanks.bam" % GENOME_FLANKS_DIR,
        genome_flanks_mate_bam = "%s/{sample}.{genome}.{insertseq}.flankmates.bam" % GENOME_FLANKS_DIR,
        genome_noflanks_bam= "%s/{sample}.{genome}.{insertseq}.noflanks.bam" % GENOME_FLANKS_DIR,
        insertseq_flanks_bam = "%s/{sample}.{genome}.{insertseq}.flanks.bam" % INSERTSEQ_FLANKS_DIR,
        insertseq_flanks_mate_bam = "%s/{sample}.{genome}.{insertseq}.flankmates.bam" % INSERTSEQ_FLANKS_DIR,
        insertseq_noflanks_bam = "%s/{sample}.{genome}.{insertseq}.noflanks.bam" % INSERTSEQ_FLANKS_DIR
    params:
        genome_flanks_sam = "%s/{sample}.{genome}.{insertseq}.flanks.sam" % GENOME_FLANKS_DIR,
        genome_flanks_mate_sam = "%s/{sample}.{genome}.{insertseq}.flankmates.sam" % GENOME_FLANKS_DIR,
        genome_noflanks_sam= "%s/{sample}.{genome}.{insertseq}.noflanks.sam" % GENOME_FLANKS_DIR,
        insertseq_flanks_sam = "%s/{sample}.{genome}.{insertseq}.flanks.sam" % INSERTSEQ_FLANKS_DIR,
        insertseq_flanks_mate_sam = "%s/{sample}.{genome}.{insertseq}.flankmates.sam" % INSERTSEQ_FLANKS_DIR,
        insertseq_noflanks_sam = "%s/{sample}.{genome}.{insertseq}.noflanks.sam" % INSERTSEQ_FLANKS_DIR,
        sample='{sample}'
    shell:
        "python scripts/get_flanking_reads.py {input.genome_sam} {input.insertseq_sam} "
        "{input.class1} {input.class2} {params.genome_flanks_sam} {params.genome_flanks_mate_sam} {params.genome_noflanks_sam} "
        "{params.insertseq_flanks_sam} {params.insertseq_flanks_mate_sam} {params.insertseq_noflanks_sam} {params.sample}; "
        "samtools sort {params.genome_flanks_sam} > {output.genome_flanks_bam}; "
        "samtools sort {params.genome_flanks_mate_sam} > {output.genome_flanks_mate_bam}; "
        "samtools sort {params.genome_noflanks_sam} > {output.genome_noflanks_bam}; "
        "samtools sort {params.insertseq_flanks_sam} > {output.insertseq_flanks_bam}; "
        "samtools sort {params.insertseq_flanks_mate_sam} > {output.insertseq_flanks_mate_bam}; "
        "samtools sort {params.insertseq_noflanks_sam} > {output.insertseq_noflanks_bam}; "
        "rm {params.genome_flanks_sam} {params.genome_flanks_mate_sam} {params.genome_noflanks_sam} "
        "{params.insertseq_flanks_sam} {params.insertseq_flanks_mate_sam} {params.insertseq_noflanks_sam};"

rule samtools_index_flank_bams:
    input:
        genome_flanks_bam = "%s/{sample}.{genome}.{insertseq}.flanks.bam" % GENOME_FLANKS_DIR,
        genome_flanks_mate_bam = "%s/{sample}.{genome}.{insertseq}.flankmates.bam" % GENOME_FLANKS_DIR,
        genome_noflanks_bam= "%s/{sample}.{genome}.{insertseq}.noflanks.bam" % GENOME_FLANKS_DIR,
        insertseq_flanks_bam = "%s/{sample}.{genome}.{insertseq}.flanks.bam" % INSERTSEQ_FLANKS_DIR,
        insertseq_flanks_mate_bam = "%s/{sample}.{genome}.{insertseq}.flankmates.bam" % INSERTSEQ_FLANKS_DIR,
        insertseq_noflanks_bam = "%s/{sample}.{genome}.{insertseq}.noflanks.bam" % INSERTSEQ_FLANKS_DIR
    output:
        genome_flanks_bam_bai = "%s/{sample}.{genome}.{insertseq}.flanks.bam.bai" % GENOME_FLANKS_DIR,
        genome_flanks_mate_bam_bai = "%s/{sample}.{genome}.{insertseq}.flankmates.bam.bai" % GENOME_FLANKS_DIR,
        genome_noflanks_bam_bai= "%s/{sample}.{genome}.{insertseq}.noflanks.bam.bai" % GENOME_FLANKS_DIR,
        insertseq_flanks_bam_bai = "%s/{sample}.{genome}.{insertseq}.flanks.bam.bai" % INSERTSEQ_FLANKS_DIR,
        insertseq_flanks_mate_bam_bai = "%s/{sample}.{genome}.{insertseq}.flankmates.bam.bai" % INSERTSEQ_FLANKS_DIR,
        insertseq_noflanks_bam_bai = "%s/{sample}.{genome}.{insertseq}.noflanks.bam.bai" % INSERTSEQ_FLANKS_DIR
    shell:
        "samtools index {input.genome_flanks_bam}; "
        "samtools index {input.genome_flanks_mate_bam}; "
        "samtools index {input.genome_noflanks_bam}; "
        "samtools index {input.insertseq_flanks_bam}; "
        "samtools index {input.insertseq_flanks_mate_bam}; "
        "samtools index {input.insertseq_noflanks_bam};"

