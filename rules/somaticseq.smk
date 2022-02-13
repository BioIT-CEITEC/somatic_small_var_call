import re

if config["is_paired"] == True:
    available_varcallers = {
        "muse":"single_output",
        "mutect2":"single_output",
        "scalpel":"single_output",
        "somaticsniper":"single_output",
        "vardict":"single_output",
        "lofreq_paired":"multi_output",
        "strelka_paired":"multi_output",
        "varscan_paired":"multi_output",
    }
else:
    available_varcallers = {
        "mutect2":"single_output",
        "scalpel":"single_output",
        "vardict":"single_output",
        "lofreq_single":"single_output",
        "strelka_single":"single_output",
        "varscan_single":"single_output",
    }

# variant callers that have separate output files for SNVs and Indels


def individual_caller_outputs(wildcards):
    output_list = [expand(eval('.'.join(['rules', varcaller, 'output'])),sample_name=wildcards.sample_name)
     for varcaller in available_varcallers.keys() if any([varcaller.find(config_caller) != -1 for config_caller in callers])]
    flat_list = [item for sublist in output_list for item in sublist]
    return(flat_list)

def inputflags(wildcards):
    somaticseq_inputflags = []
    for varcaller in available_varcallers.keys():
        varcaller_short = re.sub("_paired|_single","",varcaller)
        if varcaller_short in callers:
            if available_varcallers[varcaller] == "single_output":
                somaticseq_inputflags += ["--" + varcaller_short + "-vcf"]
            else:
                somaticseq_inputflags += ["--" + varcaller_short + "-snv"]
                somaticseq_inputflags += ["--" + varcaller_short + "-indel"]

    caller_output = individual_caller_outputs(wildcards)

    return(expand("{somaticseq_inputflag} {caller_output}",zip\
        ,somaticseq_inputflag=somaticseq_inputflags\
        ,caller_output=caller_output \
        ))

rule somaticseq:
    input:
        unpack(bam_inputs),
        caller_output = individual_caller_outputs,
        ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir = reference_directory,ref_name = config["reference"])[0],
        regions =  expand("{ref_dir}/intervals/{library_scope}/{library_scope}.bed",ref_dir = reference_directory,library_scope = config["lib_ROI"])[0],
        dbsnp = expand("{ref_dir}/annot/dbSNP/common_all.vcf.gz",ref_dir=reference_directory)[0]
    output:
        snv="somatic_seq_results/{sample_name}/Consensus.sSNV.vcf",
        indel="somatic_seq_results/{sample_name}/Consensus.sINDEL.vcf"
    log: "logs/{sample_name}/somaticseq.log"
    threads: 20
    resources:
        mem_mb=8000
    params:
        mincallers = 0.4,
        calling_type = config["is_paired"],
        inputflags = inputflags,
    conda:  "../wrappers/somaticseq/env.yaml"
    script: "../wrappers/somaticseq/script.py"


rule postprocess_somaticseq_variants:
    input:
        snv="somatic_seq_results/{sample_name}/Consensus.sSNV.vcf",
        indel="somatic_seq_results/{sample_name}/Consensus.sINDEL.vcf"
    output:
        var_tab = "somatic_seq_results/{sample_name}.variants.tsv"
    log: "logs/{sample_name}/postprocess_somaticseq_variants.log"
    threads: 1
    resources:
        mem_mb=8000
    params:
        calling_type=config["is_paired"]
    conda:  "../wrappers/postprocess_somaticseq_variants/env.yaml"
    script: "../wrappers/postprocess_somaticseq_variants/script.py"
