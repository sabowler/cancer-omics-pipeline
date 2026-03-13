"""Trimmomatic — paired-end adapter trimming and quality filtering."""

rule trimmomatic:
    input:
        r1 = str(PROJ / "fastq/{sample}_R1.fastq"),
        r2 = str(PROJ / "fastq/{sample}_R2.fastq"),
    output:
        r1_trimmed   = str(PROJ / "fastq/trimmed/{sample}_trimmed_1.fastq"),
        r2_trimmed   = str(PROJ / "fastq/trimmed/{sample}_trimmed_2.fastq"),
        r1_unpaired  = str(PROJ / "fastq/unpaired/{sample}_unpaired_1.fastq"),
        r2_unpaired  = str(PROJ / "fastq/unpaired/{sample}_unpaired_2.fastq"),
    params:
        adapters  = config["references"]["adapters"],
        leading   = config["trimmomatic"]["leading"],
        trailing  = config["trimmomatic"]["trailing"],
        avgqual   = config["trimmomatic"]["avgqual"],
        minlen    = config["trimmomatic"]["minlen"],
    threads: config["trimmomatic"]["threads"]
    log:
        str(PROJ / "logs/trimmomatic/{sample}.log")
    conda: "../envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1_trimmed} {output.r1_unpaired} \
            {output.r2_trimmed} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            AVGQUAL:{params.avgqual} \
            MINLEN:{params.minlen} \
            2> {log}
        """
