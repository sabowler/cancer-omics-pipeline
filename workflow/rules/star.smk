"""STAR — splice-aware RNA-seq alignment."""

rule star_align:
    input:
        r1 = str(PROJ / "fastq/trimmed/{sample}_trimmed_1.fastq"),
        r2 = str(PROJ / "fastq/trimmed/{sample}_trimmed_2.fastq"),
    output:
        bam = str(PROJ / "bam/{sample}.bam"),
        bai = str(PROJ / "bam/{sample}.bam.bai"),
    params:
        genome_dir = config["references"]["star_index"],
        gtf        = config["references"]["gtf"],
        prefix     = str(PROJ / "bam/star_tmp/{sample}/"),
        multimap   = config["star"]["out_filter_multimap_nmax"],
        matchmin   = config["star"]["out_filter_match_nmin"],
    threads: config["star"]["threads"]
    log:
        str(PROJ / "logs/star/{sample}.log")
    conda: "../envs/star.yaml"
    shell:
        """
        mkdir -p {params.prefix}
        STAR \
            --runThreadN {threads} \
            --runMode alignReads \
            --outFilterMultimapNmax {params.multimap} \
            --outFilterMatchNmin {params.matchmin} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --twopassMode Basic \
            --outFileNamePrefix {params.prefix} \
            --genomeDir {params.genome_dir} \
            --sjdbGTFfile {params.gtf} \
            --readFilesIn {input.r1} {input.r2} \
            2> {log}
        mv {params.prefix}*sortedByCoord.out.bam {output.bam}
        samtools index {output.bam}
        rm -rf {params.prefix}
        """
