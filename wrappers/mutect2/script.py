######################################
# wrapper for rule: mutect2
######################################
import os
import subprocess
from snakemake.shell import shell


log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: mutect2 \n##\n")
f.close()


shell.executable("/bin/bash")

# ZAKOMENTOVANO 23.11.2020 small variants analysis
# normal_cfg = snakemake.params.normal_sample
# normal_sample_name = str(normal_cfg.loc[normal_cfg.fullname == snakemake.wildcards.fullname,"sample"].min())
# if normal_sample_name == "nan":
#     normal_sample_name = normal_cfg["sample"].min()

version = str(subprocess.Popen("gatk Mutect2 --version true 2>&1 | grep \"[Vv]ersion:\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'a+')
f.write("## VERSION: gatk "+version+"\n")
f.close()

if os.stat(snakemake.input.regions).st_size == 0:
    intervals_call = ""
else:
    intervals_call = " -L " + snakemake.input.regions

command = "mkdir -p " + os.path.dirname(snakemake.params.bamout)
# command = "mkdir -p " + os.path.dirname(''.join(snakemake.params.bamout))
shell(command)

non_filtered_vcf = snakemake.output.vcf.replace(".vcf",".not_filtered.vcf")

if snakemake.params.calling_type:
    command = "gatk --java-options \"-Xmx30g \" Mutect2" + \
              " -R " + snakemake.input.ref + \
              " -I " + snakemake.input.tumor + \
              " -I " + snakemake.input.normal + \
              " -tumor " + str(snakemake.params.sample_orig_bam_names["tumor"]) + \
              " -normal " + str(snakemake.params.sample_orig_bam_names["normal"]) + \
              " --native-pair-hmm-threads " + str(snakemake.threads) + \
              intervals_call + \
              " -O " + non_filtered_vcf + \
              " -bamout " + snakemake.params.bamout + \
              " >> " + log_filename + " 2>&1"
else:
    command = "gatk --java-options \"-Xmx30g \" Mutect2" + \
              " -R " + snakemake.input.ref + \
              " -I " + snakemake.input.tumor + \
              " -tumor " + str(snakemake.params.sample_orig_bam_names["tumor"]) + \
              " --native-pair-hmm-threads " + str(snakemake.threads) + \
              intervals_call + \
              " -O " + non_filtered_vcf + \
              " -bamout " + snakemake.params.bamout + \
              " >> " + log_filename + " 2>&1"



f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

command = "gatk --java-options \"-Xmx30g\" FilterMutectCalls" + \
               " -R "+ snakemake.input.ref +\
               " -V "+ non_filtered_vcf +\
               " -O " + snakemake.output.vcf +\
               " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)
