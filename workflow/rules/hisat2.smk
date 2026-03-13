"""HISAT2 — alternative splice-aware RNA-seq aligner."""

rule hisat2_align:
    input:
        r1 = str(PROJ / "fastq/trimmed/{sample}_trimmed_1.fastq"),
        r2 = str(PROJ / "fastq/trimmed/{sample}_trimmed_2.fastq"),
    output:
        bam = str(PROJ / "bam/{sample}.bam"),
        bai = str(PROJ / "bam/{sample}.bam.bai"),
    params:
        index = config["references"]["hisat2_index"],
    threads: config["hisat2"]["threads"]
    log:
        str(PROJ / "logs/hisat2/{sample}.log")
    conda: "../envs/hisat2.yaml"
    shell:
        """
        hisat2 -p {threads} --dta \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --summary-file {log} \
            | samtools sort -o {output.bam}
        samtools index {output.bam}
        """
