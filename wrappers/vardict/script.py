######################################
# wrapper for rule: vardict
######################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

log_filename = str(snakemake.log)

if snakemake.params.calling_type:
    TEST_STRAND_BIAS = os.path.abspath(os.path.dirname(__file__))+"/testsomatic.R"
else:
    TEST_STRAND_BIAS = os.path.abspath(os.path.dirname(__file__)) + "/teststrandbias.R"

#changed in newer VARDICT versions, VAR2VCF is different for single/paired samples
if snakemake.params.calling_type:
    VAR2VCF = os.path.abspath(os.path.dirname(__file__))+"/var2vcf_paired.pl"
else:
    VAR2VCF = os.path.abspath(os.path.dirname(__file__))+"/var2vcf_valid.pl"



f = open(log_filename, 'a+')
f.write("\n##\n## RULE: vardict \n##\n")
f.close()

shell.executable("/bin/bash")

# nic nevraci...
version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'a+')
f.write("## VERSION: vardict-java 1.8.2 "+version+"\n")
f.close()

#AttributeError: 'Wildcards' object has no attribute 'full_name'
# replace 'full_name' with fullname

if snakemake.params.calling_type:
    # TEST_STRAND_BIAS = os.path.abspath(os.path.dirname(__file__)) + "/testsomatic.R"
    # VAR2VCF = os.path.abspath(os.path.dirname(__file__)) + "/var2vcf_somatic.pl"

    command = "vardict-java -G " + snakemake.input.ref + \
              " -th " + str(snakemake.threads) + \
              " -N " + snakemake.wildcards.sample_name + ".tumor" + \
              " -b '"+ snakemake.input.tumor + "|" + snakemake.input.normal + "'" + \
              " -c 1 -S 2 -E 3 -g 4 " + snakemake.input.regions + " 2>> " + log_filename + \
              " | " + TEST_STRAND_BIAS + \
              " | " + VAR2VCF + \
              " -m 7 -c 1 -N '" + snakemake.wildcards.sample_name + ".tumor|" + snakemake.wildcards.sample_name + ".normal'" + \
              " -f " + str(snakemake.params.AF_threshold) + " > " + snakemake.output.vcf + ".temp1" + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "sed -i '/[ACGT]<dup/d' " + snakemake.output.vcf + ".temp1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "bcftools view -f PASS -e 'STATUS ~ \"Germline\"' " + \
              snakemake.output.vcf + ".temp1 > " + snakemake.output.vcf + ".temp2" + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "picard -Xmx10g SortVcf" + \
              " I=" + snakemake.output.vcf + ".temp2" + \
              " O=" + snakemake.output.vcf + \
              " SD=" + snakemake.input.refdict + \
              " >> " + log_filename + " 2>&1"

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    shell("rm " + snakemake.output.vcf + ".temp1")
    shell("rm " + snakemake.output.vcf + ".temp2")
else:
    # TEST_STRAND_BIAS = os.path.abspath(os.path.dirname(__file__)) + "/teststrandbias.R"
    # VAR2VCF = os.path.abspath(os.path.dirname(__file__)) + "/var2vcf_valid.pl"

    command = "vardict-java -G " + snakemake.input.ref + \
              " -th " + str(snakemake.threads) + \
              " -N " + snakemake.wildcards.sample_name + \
              " -b " + snakemake.input.tumor + \
              " -c 1 -S 2 -E 3 -g 4 " + snakemake.input.regions + " 2>> " + log_filename + \
              " | " + TEST_STRAND_BIAS + \
              " | " + VAR2VCF + \
              " -m 7 -c 1 -N " + snakemake.wildcards.sample_name + \
              " -f " + str(snakemake.params.AF_threshold) + " > " + snakemake.output.vcf + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "sed -i '/[ACGT]<dup/d' " + snakemake.output.vcf

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

#  removed unused logs
    # command = "vardict-java -G " + snakemake.input.ref + \
    #           " -th " + str(snakemake.threads) + \
    #           " -N " + snakemake.wildcards.sample_name + \
    #           " -b " + snakemake.input.tumor + \
    #           " -c 1 -S 2 -E 3 -g 4 " + snakemake.input.regions + \
    #           " 2>> " + log_filename + \
    #           " | " + TEST_STRAND_BIAS + \
    #           " 2>> " + log_filename + \
    #           " | " + VAR2VCF + \
    #           " -m 7 -c 1 -N " + snakemake.wildcards.sample_name + \
    #           " -f " + str(snakemake.params.AF_threshold) + " > " + snakemake.output.vcf + \
    #           " 2>> " + log_filename


    # command = "vardict-java -G " + snakemake.input.ref + \
    #           " -th " + str(snakemake.threads) + \
    #           " -N " + snakemake.wildcards.sample_name + ".tumor" + \
    #           " -b '"+ snakemake.input.tumor + "|" + snakemake.input.normal + "'" + \
    #           " -c 1 -S 2 -E 3 -g 4 " + snakemake.input.regions + " 2>> " + log_filename + \
    #           " | " + TEST_STRAND_BIAS + \
    #           " | " + VAR2VCF + \
    #           " -m 7 -c 1 -N '" + snakemake.wildcards.sample_name + ".tumor|" + snakemake.wildcards.sample_name + ".normal'" + \
    #           " -f " + str(snakemake.params.AF_threshold) + " > " + snakemake.output.vcf + ".temp1" + \
    #           " 2>> " + log_filename