"""Pathoscope — microbiome profiling via competitive read mapping."""

PATHOSCOPE_REFS = config["pathoscope"]["references"]

rule pathoscope_map:
    input:
        r1 = str(PROJ / "fastq/trimmed/{sample}_trimmed_1.fastq"),
        r2 = str(PROJ / "fastq/trimmed/{sample}_trimmed_2.fastq"),
    output:
        report = str(PROJ / "Pathoscope_GRCh38.CHM13v2/{sample}-sam-report.tsv"),
    params:
        microbiome_dir  = config["pathoscope"]["microbiome_dir"],
        target_index    = config["pathoscope"]["target_index"],
        filter_prefixes = config["pathoscope"]["filter_prefixes"],
        theta_prior     = config["pathoscope"]["theta_prior"],
        outdir          = str(PROJ / "Pathoscope_GRCh38.CHM13v2"),
    threads: config["pathoscope"]["threads"]
    log:
        str(PROJ / "logs/pathoscope/{sample}.log")
    conda: "../envs/pathoscope.yaml"
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.sample}
        cd {params.microbiome_dir}

        pathoscope MAP \
            -1 {input.r1} -2 {input.r2} \
            -targetIndexPrefixes {params.target_index} \
            -filterIndexPrefixes {params.filter_prefixes} \
            -outDir {params.outdir}/{wildcards.sample} \
            -outAlign {wildcards.sample}.sam \
            -expTag {wildcards.sample} 2> {log}

        pathoscope ID \
            -alignFile {params.outdir}/{wildcards.sample}/{wildcards.sample}.sam \
            -fileType sam \
            -outDir {params.outdir}/{wildcards.sample} \
            -expTag {wildcards.sample} \
            -thetaPrior {params.theta_prior} 2>> {log}

        mv {params.outdir}/{wildcards.sample}/{wildcards.sample}-sam-report.tsv {output.report}
        rm -rf {params.outdir}/{wildcards.sample}
        """

rule merge_pathoscope:
    input:
        expand(str(PROJ / "Pathoscope_GRCh38.CHM13v2/{sample}-sam-report.tsv"), sample=ALL_SAMPLES)
    output:
        str(PROJ / "Pathoscope_GRCh38.CHM13v2/Combined.tsv")
    log:
        str(PROJ / "logs/pathoscope/merge.log")
    conda: "../envs/base_python.yaml"
    script:
        "../scripts/merge_pathoscope.py"
