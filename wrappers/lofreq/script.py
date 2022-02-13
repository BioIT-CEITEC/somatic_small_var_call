######################################
# wrapper for rule: lofreq
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## WRAPPER: lofreq \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

if snakemake.params.calling_type:
    command = "lofreq somatic " + \
              " --threads " + str(snakemake.threads) + \
              " -t " + snakemake.input.tumor + \
              " -n " + snakemake.input.normal + \
              " --call-indels " + \
              " -f " + snakemake.input.ref + \
              " -l " + snakemake.input.regions + \
              " -o " + snakemake.params.prefix + \
              " -d " + snakemake.input.dbsnp + \
              " >> " + log_filename + " 2>&1"
else:
    command = "lofreq call " + \
              " --call-indels " + \
              " -f " + snakemake.input.ref + \
              " -l " + snakemake.input.regions + \
              " -o " + snakemake.output.vcf + \
              " -d " + snakemake.input.dbsnp + \
              " " + snakemake.input.tumor + \
              " >> " + log_filename + " 2>&1"



f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
