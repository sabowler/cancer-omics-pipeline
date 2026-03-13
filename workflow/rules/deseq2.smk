"""DESeq2 — differential expression analysis from HTSeq count matrices."""

rule deseq2:
    input:
        counts  = expand(str(PROJ / "htseq/{sample}.tsv"), sample=ALL_SAMPLES),
        samples = SAMPLES_TSV,
    output:
        results = str(PROJ / "deseq2/results.xlsx"),
        rdata   = str(PROJ / "deseq2/dds.RData"),
    params:
        condition_col   = config["deseq2"]["condition_col"],
        reference_level = config["deseq2"]["reference_level"],
        htseq_dir       = str(PROJ / "htseq"),
        out_dir         = str(PROJ / "deseq2"),
    log:
        str(PROJ / "logs/deseq2/deseq2.log")
    conda: "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"
