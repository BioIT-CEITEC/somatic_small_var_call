
# ####################################
# # AFTER ANNOTATION PROCESSING
# #

rule process_and_format_annot_variants:
    input:  var_tabs = expand("somatic_seq_results/{sample_name}.variants.tsv", sample_name = sample_tab.sample_name),
            annotated = "annotate/all_variants.annotated.processed.tsv",
            format_file = workflow.basedir + "/resources/formats/" + config["format"] + ".txt",
    output: all_vars_xlsx = "final_variant_table.xlsx",
            all_vars_tsv = "final_variant_table.tsv",
            mut_loads = "mutation_loads.xlsx",
            per_sample_var_tabs = expand("per_sample_final_var_tabs/{sample_name}.variants.xlsx", sample_name = sample_tab.sample_name),
    log:    "logs/postprocess_and_format_annot_variants.log"
    threads: 10
    resources:
        mem_mb=8000
    params: reference = config["reference"],
            min_variant_frequency = str(config["min_variant_frequency"]),
            format = config["format"],
            anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"]),
            ref_dir = reference_directory
    conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
    script: "../wrappers/process_and_format_annot_variants/script.py"


# rule cohort_related_processing:
#     input:  var_tabs = expand("somatic_seq_results/{sample_name}.variants.tsv", sample_name = sample_tab.sample_name),
#             annotated = "annotate/all_variants.annotated.tsv"
#     output: all_vars = "cohort_data"
#     log:    "logs/cohort_related_processing.log"
#     threads: 10
#     params: reference = config["reference"],
#             min_variant_frequency = str(config["min_variant_frequency"]),
#             format = config["format"],
#             anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"])
#     conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
#     script: "../wrappers/process_and_format_annot_variants/script.py"
#
# rule cohort_sample_check:
#     input:  var_tabs = expand("somatic_seq_results/{sample_name}.variants.tsv", sample_name = sample_tab.sample_name),
#             annotated = "annotate/all_variants.annotated.tsv"
#     output: all_vars = "full_variant_table.xlsx"
#     log:    "logs/postprocess_and_format_annot_variants.log"
#     threads: 10
#     params: reference = config["reference"],
#             min_variant_frequency = str(config["min_variant_frequency"]),
#             format = config["format"],
#             anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"])
#     conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
#     script: "../wrappers/process_and_format_annot_variants/script.py"