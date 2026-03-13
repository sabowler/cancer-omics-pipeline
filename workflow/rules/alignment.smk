"""
Rules: alignment
RNA-seq alignment (STAR or HISAT2) and DNA alignment (Bowtie2 + GATK BQSR).
Aligner is selected via config["rna_aligner"]: "star" | "hisat2"
"""

# ── RNA: STAR ──────────────────────────────────────────────────────────────────
rule star_align:
    input:
        r1 = OUTDIR + "/trimmed/{sample}_trimmed_1.fastq",
        r2 = OUTDIR + "/trimmed/{sample}_trimmed_2.fastq",
    output:
        bam     = OUTDIR + "/bam/{sample}.bam",
        counts  = OUTDIR + "/bam/{sample}ReadsPerGene.out.tab",
    log:
        OUTDIR + "/logs/star/{sample}.log"
    conda:
        "../envs/alignment.yaml"
    params:
        star_index   = config["references"]["star_index"],
        gtf          = config["references"]["gtf"],
        threads      = config["star"]["threads"],
        multimap     = config["star"]["out_filter_multimap"],
        matchmin     = config["star"]["out_filter_match_min"],
        tmp_dir      = config["tmp_dir"] + "/star_{sample}",
        prefix       = OUTDIR + "/bam/{sample}",
    shell:
        """
        STAR \
            --runThreadN {params.threads} \
            --runMode alignReads \
            --outFilterMultimapNmax {params.multimap} \
            --outFilterMatchNmin {params.matchmin} \
            --outSAMtype BAM Unsorted SortedByCoordinate \
            --quantMode GeneCounts \
            --twopassMode Basic \
            --outFileNamePrefix {params.prefix} \
            --genomeDir {params.star_index} \
            --sjdbGTFfile {params.gtf} \
            --readFilesIn {input.r1} {input.r2} \
            2> {log}
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam}
        samtools index {output.bam}
        """


# ── RNA: HISAT2 ────────────────────────────────────────────────────────────────
rule hisat2_align:
    input:
        r1 = OUTDIR + "/trimmed/{sample}_trimmed_1.fastq",
        r2 = OUTDIR + "/trimmed/{sample}_trimmed_2.fastq",
    output:
        bam = OUTDIR + "/bam/{sample}.bam",
    log:
        summary = OUTDIR + "/logs/hisat2/{sample}.txt",
    conda:
        "../envs/alignment.yaml"
    params:
        index   = config["references"]["hisat2_index"],
        threads = config["hisat2"]["threads"],
        tmp_sam = OUTDIR + "/bam/{sample}.sam",
    shell:
        """
        hisat2 -p {params.threads} --dta \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            -S {params.tmp_sam} \
            --summary-file {log.summary}
        samtools sort -o {output.bam} {params.tmp_sam}
        samtools index {output.bam}
        rm {params.tmp_sam}
        """


# ── DNA: Bowtie2 + GATK BQSR ──────────────────────────────────────────────────
rule bowtie2_align:
    input:
        r1 = OUTDIR + "/trimmed/{sample}_trimmed_1.fastq",
        r2 = OUTDIR + "/trimmed/{sample}_trimmed_2.fastq",
    output:
        bam = OUTDIR + "/bam/{sample}.bam",
    log:
        OUTDIR + "/logs/bowtie2/{sample}.log"
    conda:
        "../envs/alignment.yaml"
    params:
        index   = config["references"]["genome_index"],
        fasta   = config["references"]["genome_fasta"],
        gnomad  = config["references"]["gnomad_vcf"],
        threads = config["bowtie2"]["threads"],
        metrics = OUTDIR + "/bam/metrics/{sample}.metrics.txt",
        table   = OUTDIR + "/bam/tables/{sample}_recal_data.table",
    shell:
        """
        bowtie2 -p {params.threads} \
            --rg-id 1 --rg LB:lib1 --rg PL:ILLUMINA --rg PU:unit1 --rg SM:{wildcards.sample} \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} | \
        samtools view -bS | \
        samtools sort -@ {params.threads} -o {OUTDIR}/bam/{wildcards.sample}.raw.bam

        gatk MarkDuplicates \
            -I {OUTDIR}/bam/{wildcards.sample}.raw.bam \
            -O {OUTDIR}/bam/{wildcards.sample}.marked.bam \
            -M {params.metrics}

        gatk BaseRecalibrator \
            -I {OUTDIR}/bam/{wildcards.sample}.marked.bam \
            -R {params.fasta} \
            --known-sites {params.gnomad} \
            -O {params.table}

        gatk ApplyBQSR \
            -R {params.fasta} \
            -I {OUTDIR}/bam/{wildcards.sample}.marked.bam \
            --bqsr-recal-file {params.table} \
            -O {output.bam}

        samtools index {output.bam}
        rm {OUTDIR}/bam/{wildcards.sample}.raw.bam {OUTDIR}/bam/{wildcards.sample}.marked.bam
        2> {log}
        """


# ── Select RNA aligner at runtime ──────────────────────────────────────────────
def get_alignment_rule(wildcards):
    """Route to STAR or HISAT2 based on config."""
    aligner = config.get("rna_aligner", "star").lower()
    if aligner == "hisat2":
        return rules.hisat2_align.output.bam.format(**wildcards)
    return rules.star_align.output.bam.format(**wildcards)
