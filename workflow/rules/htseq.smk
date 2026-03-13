"""HTSeq — gene-level read counting from BAM files."""

rule htseq_count:
    input:
        bam = str(PROJ / "bam/{sample}.bam"),
        gtf = config["references"]["gtf"],
    output:
        counts = str(PROJ / "htseq/{sample}.tsv"),
    params:
        fmt      = config["htseq"]["format"],
        order    = config["htseq"]["order"],
        mode     = config["htseq"]["mode"],
        stranded = config["htseq"]["stranded"],
        minaqual = config["htseq"]["minaqual"],
        feat     = config["htseq"]["type"],
        idattr   = config["htseq"]["idattr"],
    log:
        str(PROJ / "logs/htseq/{sample}.log")
    conda: "../envs/htseq.yaml"
    shell:
        """
        htseq-count \
            --format {params.fmt} \
            --order {params.order} \
            --mode {params.mode} \
            --stranded {params.stranded} \
            --minaqual {params.minaqual} \
            --type {params.feat} \
            --idattr {params.idattr} \
            {input.bam} {input.gtf} \
            > {output.counts} 2> {log}
        """
