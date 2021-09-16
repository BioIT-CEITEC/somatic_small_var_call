######################################
# wrapper for rule: germline_strelka
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: somatic strelka \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("configureStrelkaSomaticWorkflow.py --version 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: strelka "+version+"\n")
f.close()

#not implemented --targeted
if snakemake.params.library_scope == "wgs" or snakemake.params.sample_material == "RNA":
    scope = ""
else:
    scope = " --exome --callRegions " + snakemake.input.regions_gz + " "


shell("rm -fR " + snakemake.params.dir)

if snakemake.params.calling_type == "paired":
    command = "configureStrelkaSomaticWorkflow.py" + \
              scope + \
              " --normalBam " + snakemake.input.normal + \
              " --tumorBam " + snakemake.input.tumor + \
              " --referenceFasta " + snakemake.input.ref + \
              " --runDir " + snakemake.params.dir + \
              " --disableEVS " + \
              " >> " + log_filename + " 2>&1"
else:
    command = "configureStrelkaGermlineWorkflow.py" + \
              scope + \
              " --bam " + snakemake.input.tumor + \
              " --referenceFasta " + snakemake.input.ref + \
              " --runDir " + snakemake.params.dir + \
              " >> " + log_filename + " 2>&1"



f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = snakemake.params.dir + "/runWorkflow.py -m local -j " + str(snakemake.threads) + " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# remove workspace - lots of files messing copying speeds
shell("rm -fR " + snakemake.params.dir + "/workspace")






# command = "gatk MergeVcfs " + \
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz "+\
#         " -I "+ snakemake.params.dir + "/results/variants/somatic.snvs.vcf.gz "+\
#         " -R "+ snakemake.input.ref +\
#         " -D " + snakemake.input.dict +\
#         " -O " + snakemake.output.vcf +\
#         " >> " + log_filename + " 2>&1"
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)


# command =  "cp " + snakemake.params.dir + "/results/variants/somatic.indels.vcf.gz " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# command =  "gunzip " + str(snakemake.output.vcf) + ".gz"
#
# f = open(log_filename, 'a+')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
