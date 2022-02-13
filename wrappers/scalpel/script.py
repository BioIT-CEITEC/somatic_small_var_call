######################################
# wrapper for rule: scalpel
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

vcfsorter = os.path.abspath(os.path.dirname(__file__))+"/vcfsorter.pl"

f = open(log_filename, 'a+')
f.write("\n##\n## WRAPPER: scalpel \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

if snakemake.params.calling_type:
    command = "scalpel-discovery --somatic " + \
              " --ref " + snakemake.input.ref + \
              " --bed " + snakemake.input.regions + \
              " --tumor " + snakemake.input.tumor + \
              " --normal " + snakemake.input.normal + \
              " --window 600" + \
              " --dir " + snakemake.params.dir + \
              " >> " + log_filename + " 2>&1" + \
              " && scalpel-export --somatic " + \
              " --ref " + snakemake.input.ref + \
              " --bed " + snakemake.input.regions + \
              " --db " + snakemake.params.db + \
              " 2>> " + log_filename + \
              " | " + vcfsorter + \
              " -d " + snakemake.input.refdict + \
              " -v /dev/stdin " + \
              " -o " + snakemake.output.vcf + \
              " >> " + log_filename + " 2>&1"
else:
    command = "scalpel-discovery --single " + \
              " --ref " + snakemake.input.ref + \
              " --bed " + snakemake.input.regions + \
              " --bam " + snakemake.input.tumor + \
              " --window 600" + \
              " --dir " + snakemake.params.dir + \
              " >> " + log_filename + " 2>&1" + \
              " && scalpel-export --single " + \
              " --ref " + snakemake.input.ref + \
              " --bed " + snakemake.input.regions + \
              " --db " + snakemake.params.db + \
              " 2>> " + log_filename + \
              " | " + vcfsorter + \
              " -d " + snakemake.input.refdict + \
              " -v /dev/stdin " + \
              " -o " + snakemake.output.vcf + \
              " >> " + log_filename + " 2>&1"



f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
