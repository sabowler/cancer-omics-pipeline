"""TMB — tumor mutational burden via Mutect2, FilterMutectCalls, and Funcotator."""

rule mutect2:
    input:
        tumor  = str(PROJ / "bam/{sample}_T.bam"),
        normal = lambda wc: str(PROJ / f"bam/{get_paired_normal(wc.sample + '_T')}.bam"),
    output:
        vcf = str(PROJ / "TMB/somatic/{sample}.somatic.vcf.gz"),
    params:
        fasta  = config["references"]["genome_fasta"],
        gnomad = config["references"]["gnomad_vcf"],
        pon    = config["references"]["pon_vcf"],
    log:
        str(PROJ / "logs/tmb/{sample}.mutect2.log")
    conda: "../envs/gatk.yaml"
    shell:
        """
        gatk Mutect2 \
            -R {params.fasta} \
            -I {input.tumor} --tumor-sample {wildcards.sample}_T.bam \
            -I {input.normal} \
            --germline-resource {params.gnomad} \
            --panel-of-normals {params.pon} \
            --verbosity ERROR \
            -O {output.vcf} 2> {log}
        """

rule filter_mutect:
    input:
        vcf          = str(PROJ / "TMB/somatic/{sample}.somatic.vcf.gz"),
        tumor_pile   = str(PROJ / "TMB/tables/{sample}_T.pileups.table"),
        normal_pile  = str(PROJ / "TMB/tables/{sample}_N.pileups.table"),
    output:
        filtered_vcf = str(PROJ / "TMB/filtered/{sample}.filtered.vcf"),
        contam_table = str(PROJ / "TMB/contamination/{sample}.contamination.table"),
        segments     = str(PROJ / "TMB/segments/{sample}.segments.tsv"),
    params:
        fasta  = config["references"]["genome_fasta"],
        gnomad = config["references"]["gnomad_vcf"],
    log:
        str(PROJ / "logs/tmb/{sample}.filter.log")
    conda: "../envs/gatk.yaml"
    shell:
        """
        gatk GetPileupSummaries -R {params.fasta} \
            -I {wildcards.sample}_T.bam \
            -V {params.gnomad} -L {params.gnomad} \
            --verbosity ERROR \
            -O {input.tumor_pile}

        gatk GetPileupSummaries -R {params.fasta} \
            -I {wildcards.sample}_N.bam \
            -V {params.gnomad} -L {params.gnomad} \
            --verbosity ERROR \
            -O {input.normal_pile}

        gatk CalculateContamination \
            -I {input.tumor_pile} \
            -matched {input.normal_pile} \
            --tumor-segmentation {output.segments} \
            -O {output.contam_table}

        gatk FilterMutectCalls -R {params.fasta} \
            -V {input.vcf} \
            --contamination-table {output.contam_table} \
            --tumor-segmentation {output.segments} \
            -O {output.filtered_vcf} 2> {log}
        """

rule funcotator:
    input:
        vcf = str(PROJ / "TMB/filtered/{sample}.filtered.vcf"),
    output:
        maf = str(PROJ / "TMB/{sample}.maf"),
    params:
        fasta        = config["references"]["genome_fasta"],
        funcotator   = config["references"]["funcotator_db"],
    log:
        str(PROJ / "logs/tmb/{sample}.funcotator.log")
    conda: "../envs/gatk.yaml"
    shell:
        """
        gatk Funcotator \
            -R {params.fasta} \
            -V {input.vcf} \
            -O {output.maf} \
            --output-file-format MAF \
            --data-sources-path {params.funcotator} \
            --ref-version hg38 2> {log}
        """
