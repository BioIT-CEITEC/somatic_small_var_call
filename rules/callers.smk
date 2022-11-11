
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

# def mpileup_bam_input(wildcards):
#     if config["material"] != "RNA":
#         tag = "bam"
#     else:
#         tag = "RNAsplit.bam"
#     if wildcards.sample_pair == "tumor":
#         return expand("mapped/{tumor_bam}.{tag}",tumor_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"],tag=tag)
#     else:
#         return expand("mapped/{normal_bam}.{tag}",normal_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"],tag=tag)

def sample_orig_bam_names(wildcards):
    if config["tumor_normal_paired"] == True:
        return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"])[0], \
                'normal': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"])[0]}
    else:
        return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name"])[0]}


rule somaticsniper:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf = "variant_calls/{sample_name}/somaticsniper/SomaticSniper.vcf"
    log: "logs/{sample_name}/callers/somaticsniper.log"
    threads: 1
    resources: mem=50
    params: prefix = "variant_calls/{sample_name}/somaticsniper/"
    conda:  "../wrappers/somaticsniper/env.yaml"
    script: "../wrappers/somaticsniper/script.py"

rule lofreq_paired:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        dbsnp = expand("{ref_dir}/annot/dbSNP/common_all.vcf.gz",ref_dir=reference_directory)[0],
    output:
        snps="variant_calls/{sample_name}/lofreq/somatic_final.snvs.vcf.gz",
        indels="variant_calls/{sample_name}/lofreq/somatic_final.indels.vcf.gz"
    log: "logs/{sample_name}/callers/lofreq.log"
    threads: 5
    resources:
        mem_mb=12000
    params: prefix = "variant_calls/{sample_name}/lofreq/",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/lofreq/env.yaml"
    script: "../wrappers/lofreq/script.py"

rule lofreq_single:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        dbsnp = expand("{ref_dir}/annot/dbSNP/common_all.vcf.gz",ref_dir=reference_directory)[0],
    output:
        vcf="variant_calls/{sample_name}/lofreq/Lofreq.vcf"
    log: "logs/{sample_name}/callers/lofreq.log"
    threads: 5
    resources:
        mem_mb=12000
    params: prefix = "variant_calls/{sample_name}/lofreq/",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/lofreq/env.yaml"
    script: "../wrappers/lofreq/script.py"

rule muse:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        dbsnp = expand("{ref_dir}/annot/dbSNP/common_all.vcf.gz",ref_dir=reference_directory)[0],
    output:
        vcf = "variant_calls/{sample_name}/muse/MuSE.vcf"
    log: "logs/{sample_name}/callers/muse.log"
    threads: 1
    resources:
        mem_mb=4000
    params: dir = "variant_calls/{sample_name}/muse/MuSE",
            txt_res = "variant_calls/{sample_name}/muse/MuSE.MuSE.txt",
            threecol_bed_tmp = "variant_calls/{sample_name}/muse/threecol_regions_tmp.bed"
    conda:  "../wrappers/muse/env.yaml"
    script: "../wrappers/muse/script.py"

rule scalpel:
    input:
        unpack(bam_inputs),
        ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir = reference_directory,ref_name = config["reference"])[0],
        refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir = reference_directory,ref_name = config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf = "variant_calls/{sample_name}/scalpel/Scalpel.vcf"
    log: "logs/{sample_name}/callers/scalpel.log"
    threads: 1
    resources:
        mem_mb=32000
    params: dir = "variant_calls/{sample_name}/scalpel",
            db = "variant_calls/{sample_name}/scalpel/main/somatic.db.dir",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/scalpel/env.yaml"
    script: "../wrappers/scalpel/script.py"

rule mutect2:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf = "variant_calls/{sample_name}/mutect2/MuTect2.vcf"
    log: "logs/{sample_name}/callers/mutect2.log"
    threads: 5
    resources: mem = 30
    params: sample_orig_bam_names = sample_orig_bam_names,
            bamout = "variant_calls/{sample_name}/mutect2/realigned.bam",
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/mutect2/env.yaml"
    script: "../wrappers/mutect2/script.py"


rule strelka_paired:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions_gz=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        regions_tbi=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz.tbi",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        snps="variant_calls/{sample_name}/strelka/results/variants/somatic.snvs.vcf.gz",
        indels="variant_calls/{sample_name}/strelka/results/variants/somatic.indels.vcf.gz"
    log: "logs/{sample_name}/callers/strelka.log"
    threads: 5
    resources:
        mem_mb=6000
    params: dir = "variant_calls/{sample_name}/strelka",
            library_scope = config["lib_ROI"],
            sample_material=config["material"],
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/strelka/env.yaml"
    script: "../wrappers/strelka/script.py"

rule strelka_single:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions_gz=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
        regions_tbi=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed.gz.tbi",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf="variant_calls/{sample_name}/strelka/results/variants/variants.vcf.gz"
    log: "logs/{sample_name}/callers/strelka.log"
    threads: 5
    resources:
        mem_mb=6000
    params: dir = "variant_calls/{sample_name}/strelka",
            library_scope = config["lib_ROI"],
            sample_material=config["material"],
            calling_type = config["tumor_normal_paired"]
    conda:  "../wrappers/strelka/env.yaml"
    script: "../wrappers/strelka/script.py"

rule vardict:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf = "variant_calls/{sample_name}/vardict/VarDict.vcf"
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
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        snp="variant_calls/{sample_name}/varscan/VarScan2.snp.vcf",
        indel="variant_calls/{sample_name}/varscan/VarScan2.indel.vcf"
    log: "logs/{sample_name}/callers/varscan.log"
    threads: 1
    resources:
        mem_mb=9000
    params:
        tumor_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_tumor.mpileup.gz",
        normal_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_normal.mpileup.gz",
        extra = config["varscan_extra_params"],
        calling_type = config["tumor_normal_paired"]
    conda: "../wrappers/varscan/env.yaml"
    script: "../wrappers/varscan/script.py"

rule varscan_single:
    input:
        unpack(bam_inputs),
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
        regions=expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir=reference_directory,library_scope=config["lib_ROI"])[0],
    output:
        vcf="variant_calls/{sample_name}/varscan/VarScan2.vcf",
    log: "logs/{sample_name}/callers/varscan.log"
    threads: 1
    resources:
        mem_mb=9000
    params:
        tumor_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_tumor.mpileup.gz",
        normal_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_normal.mpileup.gz",
        snp="variant_calls/{sample_name}/varscan/VarScan2.snp.vcf",
        indel="variant_calls/{sample_name}/varscan/VarScan2.indel.vcf",
        extra = config["varscan_extra_params"],
        # " --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005",
        calling_type = config["tumor_normal_paired"]
    conda: "../wrappers/varscan/env.yaml"
    script: "../wrappers/varscan/script.py"
