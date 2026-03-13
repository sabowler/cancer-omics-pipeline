"""MSI — microsatellite instability scoring with msisensor-pro (tumor/normal pairs)."""

rule msi:
    input:
        tumor  = str(PROJ / "bam/{sample}_T.bam"),
        normal = lambda wc: str(PROJ / f"bam/{get_paired_normal(wc.sample + '_T')}.bam"),
    output:
        report = str(PROJ / "MSI/{sample}"),
    params:
        sites  = config["references"]["msi_sites"],
        outdir = str(PROJ / "MSI"),
    threads: config["msi"]["threads"]
    log:
        str(PROJ / "logs/msi/{sample}.log")
    conda: "../envs/msi.yaml"
    shell:
        """
        msisensor-pro msi \
            -d {params.sites} \
            -n {input.normal} \
            -t {input.tumor} \
            -o {output.report} \
            2> {log}
        """

rule merge_msi:
    input:
        expand(str(PROJ / "MSI/{sample}"), sample=TUMOR_SAMPLES)
    output:
        str(PROJ / "MSI/MSI.Merged.xlsx")
    log:
        str(PROJ / "logs/msi/merge.log")
    conda: "../envs/base_python.yaml"
    script:
        "../scripts/merge_msi.py"
