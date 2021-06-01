######################################
# wrapper for rule: somaticseq
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## WRAPPER: somaticseq \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

if snakemake.params.calling_type == "paired":
    call_type = "paired"
    bam_inputs = " --tumor-bam-file " + snakemake.input.tumor + " --normal-bam-file " + snakemake.input.normal
else:
    call_type = "single"
    bam_inputs = " --bam-file " + snakemake.input.tumor

command = "somaticseq_parallel.py " + \
            " --threads " + str(snakemake.threads) +\
            " --output-directory " + os.path.dirname(snakemake.output.snv) +\
            " --genome-reference " + snakemake.input.ref + \
            " --inclusion-region " + snakemake.input.regions + \
            " --minimum-num-callers " + str(snakemake.params.mincallers) +\
            " " + call_type +\
            bam_inputs +\
            " " + " ".join(snakemake.params.inputflags) +\
            " >> " + log_filename + " 2>&1"

# " --dbsnp-vcf " + snakemake.input.dbsnp + \
# " --algorithm xgboost " + \
# " --classifier-snv " + snakemake.input.classifier_snv + \
# " --classifier-indel " + snakemake.input.classifier_indel + \

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
