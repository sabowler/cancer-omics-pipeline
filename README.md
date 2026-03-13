# cancer-omics-pipeline

**End-to-end Snakemake pipeline for pan-cancer multi-omic analysis.**

Built and maintained by [Scott A. Bowler](https://github.com/sabowler) — Ndhlovu Lab, Weill Cornell Medicine.

Processes paired-end FASTQ files through QC, alignment, expression quantification, differential expression, variant calling, MSI/TMB scoring, and microbiome profiling — all parallelized across HPC or cloud infrastructure via Snakemake.

---

## Pipeline Overview

```
FASTQ → Trimmomatic → STAR / HISAT2 → HTSeq → DESeq2
                   ↘ Bowtie2 + GATK → MSI (msisensor-pro)
                                    → TMB (Mutect2 + Funcotator)
                   ↘ Pathoscope (microbiome profiling)
```

| Module | Tool | Output |
|---|---|---|
| QC trimming | Trimmomatic 0.39 | Trimmed FASTQ |
| RNA-seq alignment | STAR 2.7 / HISAT2 2.2 | BAM |
| Read counting | HTSeq 2.0 | Count TSV per sample |
| Differential expression | DESeq2 | results.xlsx |
| DNA alignment + BQSR | Bowtie2 + GATK4 | Recalibrated BAM |
| Microsatellite instability | msisensor-pro | MSI.Merged.xlsx |
| Tumor mutational burden | Mutect2 + Funcotator | MAF per sample |
| Microbiome profiling | Pathoscope 2.0 | Combined.tsv |

---

## Quickstart

```bash
git clone https://github.com/sabowler/cancer-omics-pipeline.git
cd cancer-omics-pipeline

# 1. Edit configs/config.yaml — set your project_dir and reference paths
# 2. Edit configs/samples.tsv — add your sample names and FASTQ paths

# Dry run (check workflow without executing)
snakemake --configfile configs/config.yaml --cores 16 -n

# Full run with conda environments
snakemake --configfile configs/config.yaml --cores 16 --use-conda

# HPC run via SLURM
snakemake --configfile configs/config.yaml \
    --executor slurm \
    --jobs 50 \
    --use-conda
```

---

## Configuration

Edit `configs/config.yaml` before running. Key sections:

```yaml
study_id: "phs000980"
project_dir: "/data/projects/Cancer/phs000980"

references:
  genome_fasta: "/data/references/Human/GRCh38.fa"
  star_index:   "/data/references/hg38+HIV1+GFP/STAR"
  gtf:          "/data/references/hg38+HIV1+GFP/hg38HIV1GFP.gtf"
  # ... (see configs/config.yaml for full reference list)

modules:
  trimmomatic: true
  star:        true
  hisat2:      false   # Use STAR or HISAT2, not both
  htseq:       true
  deseq2:      true
  bowtie2:     true
  msi:         true
  tmb:         true
  pathoscope:  true
```

### Sample Sheet

Edit `configs/samples.tsv`:

```
sample       type    fastq_r1                 fastq_r2                 condition
SAMPLE_001   tumor   SAMPLE_001_R1.fastq      SAMPLE_001_R2.fastq      tumor
SAMPLE_001N  normal  SAMPLE_001N_R1.fastq     SAMPLE_001N_R2.fastq     normal
```

Tumor/normal pairing is resolved automatically by sample name convention (`SAMPLE_001` → `SAMPLE_001N`).

---

## Output Structure

```
{project_dir}/
├── fastq/
│   ├── trimmed/       # Trimmomatic output
│   └── unpaired/
├── bam/               # Aligned BAM files (STAR / HISAT2 / Bowtie2+GATK)
├── htseq/             # Per-sample count files
├── deseq2/
│   ├── results.xlsx   # DESeq2 differential expression results
│   └── dds.RData      # DESeq2 dataset object
├── MSI/
│   └── MSI.Merged.xlsx
├── TMB/
│   ├── somatic/       # Raw Mutect2 VCFs
│   ├── filtered/      # FilterMutectCalls output
│   └── *.maf          # Funcotator MAF files
├── Pathoscope_GRCh38.CHM13v2/
│   ├── *-sam-report.tsv
│   └── Combined.tsv
└── logs/              # Per-rule logs
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io) ≥ 7.0
- [conda](https://docs.conda.io) or [mamba](https://github.com/mamba-org/mamba) (for `--use-conda`)
- All tool environments are defined in `workflow/envs/` and installed automatically with `--use-conda`

---

## Publications

Developed in support of PITCHER: Predictor of ImmunoTherapy response by Classifying Herv-k ExpRession (Provisional Patent #: D-11459)
Co-Inventor: Lishomwa C. Ndhlovu - Weill Cornell Medicine

## License

MIT
