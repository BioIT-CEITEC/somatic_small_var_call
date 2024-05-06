
def bam_inputs(wildcards):
    if config["material"] != "RNA":
        tag = "bam"
    else:
        tag = "RNAsplit.bam"

    if config["tumor_normal_paired"] == True:
        return {'tumor': expand("mapped/{tumor_bam}.{tag}",tumor_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"],tag=tag)[0], \
                'normal': expand("mapped/{normal_bam}.{tag}",normal_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"],tag=tag)[0]}
    else:
        return {'tumor': expand("mapped/{tumor_bam}.{tag}",tumor_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"],tag=tag)[0]}

def sample_orig_bam_names(wildcards):
    if config["tumor_normal_paired"] == True:
        return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"])[0], \
                'normal': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"])[0]}
    else:
        return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])[0]}


rule somaticsniper:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions= config["organism_dna_panel"],
    output:
        vcf = "somatic_varcalls/{sample_name}/somaticsniper/SomaticSniper.vcf"
    log: "logs/{sample_name}/callers/somaticsniper.log"
    threads: 1
    resources:
        mem_mb=6000
    params: prefix = "somatic_varcalls/{sample_name}/somaticsniper/"
    conda:  "../wrappers/somaticsniper/env.yaml"
    script: "../wrappers/somaticsniper/script.py"

rule lofreq_paired:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions= config["organism_dna_panel"],
        dbsnp = config["organism_dbsnp"],
    output:
        snps="somatic_varcalls/{sample_name}/lofreq/somatic_final.snvs.vcf.gz",
        indels="somatic_varcalls/{sample_name}/lofreq/somatic_final.indels.vcf.gz"
    log: "logs/{sample_name}/callers/lofreq.log"
    threads: 5
    resources:
        mem_mb=12000
    params: prefix = "somatic_varcalls/{sample_name}/lofreq/",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/lofreq/env.yaml"
    script: "../wrappers/lofreq/script.py"

rule lofreq_single:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions= config["organism_dna_panel"],
        dbsnp = config["organism_dbsnp"],
    output:
        vcf="somatic_varcalls/{sample_name}/lofreq/Lofreq.vcf"
    log: "logs/{sample_name}/callers/lofreq.log"
    threads: 5
    resources:
        mem_mb=12000
    params: prefix = "somatic_varcalls/{sample_name}/lofreq/",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/lofreq/env.yaml"
    script: "../wrappers/lofreq/script.py"

rule muse:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions = config["organism_dna_panel"],
        dbsnp = config["organism_dbsnp"],
    output:
        vcf = "somatic_varcalls/{sample_name}/muse/MuSE.vcf"
    log: "logs/{sample_name}/callers/muse.log"
    threads: 1
    resources:
        mem_mb=4000
    params: dir = "somatic_varcalls/{sample_name}/muse/MuSE",
            txt_res = "somatic_varcalls/{sample_name}/muse/MuSE.MuSE.txt",
            threecol_bed_tmp = "somatic_varcalls/{sample_name}/muse/threecol_regions_tmp.bed"
    conda:  "../wrappers/muse/env.yaml"
    script: "../wrappers/muse/script.py"

rule scalpel:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        refdict = config["organism_dict"],
        regions = config["organism_dna_panel"],
    output:
        vcf = "somatic_varcalls/{sample_name}/scalpel/Scalpel.vcf"
    log: "logs/{sample_name}/callers/scalpel.log"
    threads: 1
    resources:
        mem_mb=32000
    params: dir = "somatic_varcalls/{sample_name}/scalpel",
            db = "somatic_varcalls/{sample_name}/scalpel/main/somatic.db.dir",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/scalpel/env.yaml"
    script: "../wrappers/scalpel/script.py"

rule mutect2:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions = config["organism_dna_panel"],
    output:
        vcf = "somatic_varcalls/{sample_name}/mutect2/MuTect2.vcf"
    log: "logs/{sample_name}/callers/mutect2.log"
    threads: 5
    resources:
        mem_mb=6000
    params: sample_orig_bam_names = sample_orig_bam_names,
            bamout = "somatic_varcalls/{sample_name}/mutect2/realigned.bam",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/mutect2/env.yaml"
    script: "../wrappers/mutect2/script.py"


rule strelka_paired:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions_gz = config["organism_dna_panel"] + ".gz",
        regions_tbi = config["organism_dna_panel"] + ".gz.tbi",
    output:
        snps="somatic_varcalls/{sample_name}/strelka/results/variants/somatic.snvs.vcf.gz",
        indels="somatic_varcalls/{sample_name}/strelka/results/variants/somatic.indels.vcf.gz"
    log: "logs/{sample_name}/callers/strelka.log"
    threads: 5
    resources:
        mem_mb=6000
    params: dir = os.path.join(GLOBAL_TMPD_PATH,"somatic_varcalls/{sample_name}/strelka"),
            library_scope = config["lib_ROI"],
            sample_material=config["material"],
            calling_type = config["tumor_normal_paired"],
            snps = os.path.join(GLOBAL_TMPD_PATH,"somatic_varcalls/{sample_name}/strelka/results/variants/somatic.snvs.vcf.gz"),
            indels = os.path.join(GLOBAL_TMPD_PATH,"somatic_varcalls/{sample_name}/strelka/results/variants/somatic.indels.vcf.gz")
    conda:  "../wrappers/strelka/env.yaml"
    script: "../wrappers/strelka/script.py"

rule strelka_single:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions_gz = config["organism_dna_panel"] + ".gz",
        regions_tbi = config["organism_dna_panel"] + ".gz.tbi",
    output:
        vcf="somatic_varcalls/{sample_name}/strelka/results/variants/variants.vcf.gz"
    log: "logs/{sample_name}/callers/strelka.log"
    threads: 5
    resources:
        mem_mb=6000
    params: dir = os.path.join(GLOBAL_TMPD_PATH,"somatic_varcalls/{sample_name}/strelka"),
            library_scope = config["lib_ROI"],
            sample_material=config["material"],
            calling_type = config["tumor_normal_paired"],
            vcf = os.path.join(GLOBAL_TMPD_PATH,"somatic_varcalls/{sample_name}/strelka/results/variants/variants.vcf.gz"),
    conda:  "../wrappers/strelka/env.yaml"
    script: "../wrappers/strelka/script.py"

rule vardict:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        refdict = config["organism_dict"],
        regions = config["organism_dna_panel"],
    output:
        vcf = "somatic_varcalls/{sample_name}/vardict/VarDict.vcf"
    log: "logs/{sample_name}/callers/vardict.log"
    threads: 5
    resources:
        mem_mb=8000
    params:
        AF_threshold = config["min_variant_frequency"],
        calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"


rule varscan_paired:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions = config["organism_dna_panel"],
    output:
        snp="somatic_varcalls/{sample_name}/varscan/VarScan2.snp.vcf",
        indel="somatic_varcalls/{sample_name}/varscan/VarScan2.indel.vcf"
    log: "logs/{sample_name}/callers/varscan.log"
    threads: 1
    resources:
        mem_mb=9000
    params:
        tumor_pileup = "somatic_varcalls/{sample_name}/varscan/{sample_name}_tumor.mpileup.gz",
        normal_pileup = "somatic_varcalls/{sample_name}/varscan/{sample_name}_normal.mpileup.gz",
        extra = config["varscan_extra_params"],
        calling_type = config["tumor_normal_paired"]
    conda: "../wrappers/varscan/env.yaml"
    script: "../wrappers/varscan/script.py"

rule varscan_single:
    input:
        unpack(bam_inputs),
        ref = config["organism_fasta"],
        regions = config["organism_dna_panel"],
    output:
        vcf="somatic_varcalls/{sample_name}/varscan/VarScan2.vcf",
    log: "logs/{sample_name}/callers/varscan.log"
    threads: 1
    resources:
        mem_mb=9000
    params:
        tumor_pileup = "somatic_varcalls/{sample_name}/varscan/{sample_name}_tumor.mpileup.gz",
        normal_pileup = "somatic_varcalls/{sample_name}/varscan/{sample_name}_normal.mpileup.gz",
        snp="somatic_varcalls/{sample_name}/varscan/VarScan2.snp.vcf",
        indel="somatic_varcalls/{sample_name}/varscan/VarScan2.indel.vcf",
        extra = config["varscan_extra_params"],
        # " --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005",
        calling_type = config["tumor_normal_paired"]
    conda: "../wrappers/varscan/env.yaml"
    script: "../wrappers/varscan/script.py"


rule RNA_SplitNCigars:
    input: bam = "mapped/{sample_name}.bam",
           ref = config["organism_fasta"]
    output: bam = "mapped/{sample_name}.RNAsplit.bam",
            bai = "mapped/{sample_name}.RNAsplit.bam.bai",
    log:    run = "logs/{sample_name}/callers/RNA_SplitNCigars.log",
    params: bai = "mapped/{sample_name}.RNAsplit.bai"
    conda:  "../wrappers/RNA_SplitNCigars/env.yaml"
    script: "../wrappers/RNA_SplitNCigars/script.py"