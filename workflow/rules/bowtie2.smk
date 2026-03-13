"""Bowtie2 + GATK — DNA alignment with duplicate marking and base quality score recalibration."""

rule bowtie2_align:
    input:
        r1 = str(PROJ / "fastq/trimmed/{sample}_trimmed_1.fastq"),
        r2 = str(PROJ / "fastq/trimmed/{sample}_trimmed_2.fastq"),
    output:
        bam = str(PROJ / "bam/{sample}.bam"),
        bai = str(PROJ / "bam/{sample}.bam.bai"),
    params:
        index  = config["references"]["genome_index"],
        fasta  = config["references"]["genome_fasta"],
        gnomad = config["references"]["gnomad_vcf"],
        tmp    = config["tmp_dir"],
    threads: config["bowtie2"]["threads"]
    log:
        str(PROJ / "logs/bowtie2/{sample}.log")
    conda: "../envs/bowtie2.yaml"
    shell:
        """
        # Align
        bowtie2 -p {threads} \
            --rg-id 1 --rg LB:lib1 --rg PL:ILLUMINA --rg PU:unit1 --rg SM:{wildcards.sample} \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            2> {log} \
            | samtools view -bS \
            | samtools sort -@ {threads} -o {PROJ}/bam/{wildcards.sample}.raw.bam

        # Mark duplicates
        gatk MarkDuplicates \
            -I {PROJ}/bam/{wildcards.sample}.raw.bam \
            -O {PROJ}/bam/{wildcards.sample}.marked.bam \
            -M {PROJ}/bam/metrics/{wildcards.sample}.metrics.txt

        # Base quality score recalibration
        gatk BaseRecalibrator \
            -I {PROJ}/bam/{wildcards.sample}.marked.bam \
            -R {params.fasta} \
            --known-sites {params.gnomad} \
            -O {PROJ}/bam/tables/{wildcards.sample}_recal.table

        gatk ApplyBQSR \
            -R {params.fasta} \
            -I {PROJ}/bam/{wildcards.sample}.marked.bam \
            --bqsr-recal-file {PROJ}/bam/tables/{wildcards.sample}_recal.table \
            -O {output.bam}

        samtools index {output.bam}

        # Cleanup intermediates
        rm -f {PROJ}/bam/{wildcards.sample}.raw.bam \
              {PROJ}/bam/{wildcards.sample}.marked.bam
        """
