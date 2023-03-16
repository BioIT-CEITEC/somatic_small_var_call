#############################################
## ONLY FOR RNA SINGLE-END SOMATIC VARCALL ##
#############################################
# if not config["tumor_normal_paired"]:
#     available_varcallers = {
#         "mutect2":"single_output",
#         "vardict":"single_output",
#         "strelka_single":"single_output",
#         "varscan_single":"single_output",
#     }

# variant callers that have separate output files for SNVs and Indels


# def individual_caller_outputs(wildcards):
#     output_list = [expand(eval('.'.join(['rules', varcaller, 'output'])),sample_name=wildcards.sample_name)
#      for varcaller in available_varcallers.keys() if any([varcaller.find(config_caller) != -1 for config_caller in callers])]
#     flat_list = [item for sublist in output_list for item in sublist]
#     print(flat_list)
#     return(flat_list)

def individual_caller_input(wildcards):
    if wildcards.variant_caller == "strelka":
        return "somatic_varcalls/{sample_name}/strelka/results/variants/variants.vcf.gz"
    if wildcards.variant_caller == "mutect2":
        return "somatic_varcalls/{sample_name}/{variant_caller}/MuTect2.vcf"
    if wildcards.variant_caller == "vardict":
        return "somatic_varcalls/{sample_name}/{variant_caller}/VarDict.vcf"
    if wildcards.variant_caller == "varscan":
        return "somatic_varcalls/{sample_name}/{variant_caller}/VarScan2.vcf"


rule normalize_variants:
    input:  vcf = individual_caller_input,
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            dict= expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0]
    output: vcf = "somatic_varcalls/{sample_name}/{variant_caller}/{variant_caller}.norm.vcf"
    log:    "logs/{sample_name}/callers/{variant_caller}_normalization.log"
    threads: 1
    conda: "../wrappers/normalize_variants/env.yaml"
    script: "../wrappers/normalize_variants/script.py"

rule merge_variant_callers:
    input:  vcfs = lambda wildcards: expand("somatic_varcalls/{sample_name}/{variant_caller}/{variant_caller}.norm.vcf",\
                                            sample_name=wildcards.sample_name,\
                                            variant_caller = callers),
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            dict= expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: not_filtered_vcf = "somatic_varcalls/{sample_name}/merged/raw_calls.vcf",
            vcf= "somatic_varcalls/{sample_name}.final_variants.vcf",
            tsv = "somatic_varcalls/{sample_name}.final_variants.tsv"
    log:    "logs/{sample_name}/merge_variant_callers.log"
    params: min_var_reads_threshold = 2,
            min_callers_threshold = 1,
            min_variant_frequency = 0,
        # min_var_reads_threshold= config["min_var_reads_threshold"],
        # min_callers_threshold=config["min_callers_threshold"],
        # min_variant_frequency=config["min_variant_frequency"],
            tmp_dir = GLOBAL_TMPD_PATH
    threads: 1
    conda:  "../wrappers/merge_variant_callers/env.yaml"
    script: "../wrappers/merge_variant_callers/script.py"
