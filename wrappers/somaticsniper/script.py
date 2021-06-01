######################################
# wrapper for rule: somatic-sniper
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## WRAPPER: somatic-sniper \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

command = "bam-somaticsniper " + \
            " -Q 15 -s 1e-05 -F vcf " +\
            " -f "+ snakemake.input.ref +\
            " " + snakemake.input.tumor + \
            " " + snakemake.input.normal + \
            " " + snakemake.output.vcf +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
