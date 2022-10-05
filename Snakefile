import os
import pandas as pd

configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references/"
# JSON validation, only the description and parameters part, not the samples part
# from snakemake.utils import validate
#validate(config, "config.schema.json")


##### Config processing #####
#conversion from new json
if config["tumor_normal_paired"]:
    sample_tab_initial = pd.DataFrame.from_dict(config["samples"],orient="index")
    sample_tab = pd.DataFrame({"sample_name" : [],"sample_name_normal" : [],"sample_name_tumor" : []})
    sample_tab["sample_name"]=sample_tab_initial["donor"].unique()
    for index, row in sample_tab.iterrows():
        sample_tab.loc[index,"sample_name_normal"] = sample_tab_initial.loc[(sample_tab_initial["donor"]==row["sample_name"]) & (sample_tab_initial["tumor_normal"]=="normal"),"sample_name"].to_string(index=False)
        sample_tab.loc[index,"sample_name_tumor"] = sample_tab_initial.loc[(sample_tab_initial["donor"]==row["sample_name"]) & (sample_tab_initial["tumor_normal"]=="tumor"),"sample_name"].to_string(index=False)
else:
    sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")


##### Reference processing
#
config["material"] = "DNA"
if config["lib_ROI"] != "wgs" or config["lib_ROI"] != "RNA":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
else:
    if config["lib_ROI"] != "RNA":
        config["material"] = "RNA" # ??? nemelo by DNA
        config["lib_ROI"] = "wgs"

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

#### SPLIT CALLERS STRING
callers = config["callers"].split(';')


#### FOLDERS
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

####################################
# DEFAULT VALUES
if not "format" in config:
    config["format"] = "default"
if not "min_variant_frequency" in config:
    config["min_variant_frequency"] = 0
if not "not_use_merged" in config:
    config["not_use_merged"] = False

wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(sample_tab.sample_name),
    #sample_pair = "tumor|normal"

## PRINT to log the values


####################################
# SEPARATE RULES
include: "rules/callers.smk"
include: "rules/somaticseq.smk"
include: "rules/annotate.smk"
include: "rules/variant_postprocessing.smk"

####################################
# RULE ALL
somatic_only = list()

if config["tumor_normal_paired"]:
    somatic_only.append("somatic_only_final_variant_table.xlsx")


rule all:
    input:
        somatic_only,
        all_vars_xlsx = "final_variant_table.xlsx",
        all_vars_tsv = "final_variant_table.tsv",
        per_sample_var_tabs = expand("per_sample_final_var_tabs/{sample_name}.variants.xlsx", sample_name = sample_tab.sample_name)

