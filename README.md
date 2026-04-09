# BioIT-CEITEC/somatic_small_var_call

This repository provides a Snakemake-based workflow for somatic small variant calling. The pipeline integrates multiple widely used somatic variant callers and combines their results. The workflow supports both tumor–normal paired and tumor-only analyses and can be applied to DNA-seq as well as RNA-seq data. For RNA samples, spliced alignments are automatically processed using `SplitNCigarReads` prior to variant calling. The selection of variant callers and analysis mode is controlled via a central configuration file. The design allows individual variant callers to be enabled or disabled without modifying the workflow itself.

## Requirements
- Linux environment
- Snakemake ≥ 5.18.0
- Conda / Mamba
- Python

All remaining dependencies are handled by Snakemake using Conda environments defined per rule and they differ based on the selection of specific parameters and used tools.
## Parameters
**Note:** This workflow is primarily designed to be configured via an internal GUI. Most parameters have sensible defaults and should not require manual specification unless customization is needed.

### Truly Required parameters
These parameters must always be specified in *config.json*:

- `reference`  
  Reference genome used for the analysis (e.g., `GRCh38`).

- `lib_ROI`  
  Type of input data: `wgs` (whole-genome sequencing) or `rna` (RNA-seq). For RNA samples, the workflow automatically applies `SplitNCigarReads` processing.

- `globalResources`  
  Path to global resources directory containing reference genomes and metadata.

- `globalTmpdPath`  
  Path to temporary directory for intermediate files.

- `samples`  
  Dictionary mapping sample identifiers to sample metadata. Each sample must have a `sample_name` field. For tumor-normal pairs, also include `donor` (patient ID) and `tumor_normal` ("tumor" or "normal").

### Optional Metadata parameters
- `task_name`  
  Descriptive name for the analysis task.

- `entity_name`  
  Entity/project identification.

### Analysis Mode
- `tumor_normal_paired` (default: `true`)  
  Set to `true` for tumor-normal paired samples, `false` for tumor-only analysis.

### Variant Caller Configuration
Enable or disable individual callers. All default to `true`; set to `false` to skip:

- `somatic_use_strelka` — Strelka caller
- `somatic_use_vardict` — VarDict caller
- `somatic_use_mutect2` — GATK MuTect2 caller
- `somatic_use_lofreq` — LoFreq caller
- `somatic_use_varscan` — VarScan2 caller
- `somatic_use_muse` — MuSE caller (requires tumor-normal pairs)
- `somatic_use_somaticsniper` — SomaticSniper caller (requires tumor-normal pairs)

### Fine-tuning Parameters
- `varscan_extra_params` (default: `--strand-filter 0 --p-value 0.95 --min-var-freq 0.05`)  
  Additional command-line arguments for VarScan2.

- `min_variant_frequency` (default: `0`)  
  Minimum allele frequency threshold for VarDict and final variant filtering.
## Usage
The workflow is executed using Snakemake and requires a prepared configuration file and aligned BAM files for each sample. From the root directory of the repository, run:

```bash
snakemake --use-conda --cores <N>
```

### Required inputs
- `mapped/{sample}.bam`  
  BAM files containing aligned reads produced by an upstream pipeline (DNA-seq or RNA-seq).

For DNA samples, this BAM file is used directly for variant calling.

For RNA samples (`lib_ROI: rna` in the config), the workflow automatically   generates `mapped/{sample}.RNAsplit.bam` using `SplitNCigarReads`, and this file is used for all downstream variant calling steps.  
  
## Output
### Main outputs
- `somatic_varcalls/{sample_name}/`  
  Directory containing the final results for each sample.

- `somatic_varcalls/{sample_name}.final_variants.tsv`  
  Final consensus variant call set (combined from all enabled callers).

- `somatic_varcalls/{sample_name}/{caller}/`  
  Individual caller output directories (e.g., `varscan/`, `mutect2/`, `strelka/`, etc.), containing raw VCF files from each variant caller.

### Additional outputs

- `somatic_varcalls/{sample_name}.RNAsplit.bam` (RNA-seq only)  
  BAM file after SplitNCigarReads processing, used for variant calling on RNA samples.

- `logs/{sample_name}/callers/`  
  Log files from individual variant caller runs.

- `config.json` (snapshot)  
  Copy of the configuration file used for this run, stored for reproducibility.

## Repository structure
```
.
├── Snakefile                     
├── README.md                        # This file
├── rules/                          
│   ├── callers.smk                  # Rules for individual variant callers
│   └── somaticseq.smk               # Rules for SomaticSeq consensus calling
├── wrappers/                       
│   ├── RNA_SplitNCigars/     
│   │   ├── script.py
│   │   └── env.yaml
│   ├── lofreq/     
│   │   ├── script.py
│   │   └── env.yaml
│   ├── muse/         
│   │   ├── script.py
│   │   └── env.yaml
│   ├── mutect2/        
│   │   ├── script.py
│   │   └── env.yaml
│   ├── postprocess_somaticseq_variants/
│   │   ├── postprocess_somaticseq_variants.R     
│   │   ├── script.py
│   │   └── env.yaml
│   ├── scalpel/
│   │   ├── vcfsorter.pl     
│   │   ├── script.py
│   │   └── env.yaml
│   ├── somaticseq/          # Note: exact name for SomaticSeq wrapper
│   │   ├── script.py
│   │   └── env.yaml
│   ├── somaticsniper/        
│   │   ├── script.py
│   │   └── env.yaml
│   ├── strelka/        
│   │   ├── script.py
│   │   └── env.yaml
│   ├── vardict/
│   │   ├── testsomatic.R
│   │   ├── teststrandbias.R       
│   │   ├── var2vcf_paired.pl
│   │   ├── var2vcf_somatic.pl
│   │   ├── var2vcf_valid.pl
│   │   ├── script.py
│   │   └── env.yaml
│   └── varscan/
│       ├── combine_vcfs.R        
│       ├── script.py
│       └── env.yaml
└── workflow.config.json
```

