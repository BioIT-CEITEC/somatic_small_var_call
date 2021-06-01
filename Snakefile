import os
import pandas as pd

GLOBAL_REF_PATH = "/mnt/ssd/ssd_3/references"
# JSON validation, only the description and parameters part, not the samples part
# from snakemake.utils import validate
#validate(config, "config.schema.json")

####################################
# FOLDERS
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

##### Config processing #####

sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


callers = config["callers"].split(';')


# DEFAULT VALUES
if not "format_name" in config:
    config["format_name"] = "default"
if not "min_variant_frequency" in config:
    config["min_variant_frequency"] = 0
if not "not_use_merged" in config:
    config["not_use_merged"] = False

wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(sample_tab.sample_name),
    sample_pair = "tumor|normal"




####################################
# SEPARATE RULES
include: "rules/callers.smk"
include: "rules/somaticseq.smk"
include: "rules/annotate.smk"
include: "rules/variant_postprocessing.smk"

####################################
# RULE ALL
rule all:
    input:  
        all_vars_xlsx = "final_variant_table.xlsx",
        all_vars_tsv = "final_variant_table.tsv",
        per_sample_var_tabs = expand("per_sample_final_var_tabs/{sample_name}.variants.xlsx", sample_name = sample_tab.sample_name),

