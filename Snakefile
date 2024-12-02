import os
import pandas as pd

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

config = BR.load_organism()


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
    sample_tab = BR.load_sample()


##### Reference processing
#

# config = BR.load_organism()

if "material" not in config:
    config["material"] = "DNA"

# ####################################
# # create caller list from table
callers = []
if config["somatic_use_strelka"]:
    callers.append("strelka")
if config["somatic_use_vardict"]:
    callers.append("vardict")
if config["somatic_use_mutect2"]:
    callers.append("mutect2")
if config["somatic_use_lofreq"]:
    callers.append("lofreq")
if config["somatic_use_varscan"]:
    callers.append("varscan")
if config["somatic_use_muse"]:
    callers.append("muse")
if config["somatic_use_somaticsniper"]:
    callers.append("somaticsniper")

####################################
# DEFAULT VALUES
if not "min_variant_frequency" in config:
    config["min_variant_frequency"] = 0

wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(sample_tab.sample_name),
    #sample_pair = "tumor|normal"


####################################
# SEPARATE RULES
include: "rules/callers.smk"
include: "rules/somaticseq.smk"

####################################
# RULE ALL
rule all:
    input:  
       final_variants = expand("somatic_varcalls/{sample_name}.final_variants.tsv", sample_name = sample_tab.sample_name)

##### BioRoot utilities - prepare reference #####
#module PR:
#    snakefile: gitlab("bioroots/bioroots_utilities", path="prepare_reference.smk",branch="master")
#    config: config

#use rule * from PR as other_*
