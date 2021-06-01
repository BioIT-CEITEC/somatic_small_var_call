"""Snakemake wrapper for varscan somatic"""
import os
from snakemake.shell import shell
from snakemake.utils import makedirs

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: varscan \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

extra_params = snakemake.params.get("extra", "")

# Building output dirs
if snakemake.params.calling_type == "paired":
    makedirs(os.path.dirname(snakemake.output.snp))
    makedirs(os.path.dirname(snakemake.output.indel))
else:
    makedirs(os.path.dirname(snakemake.output.vcf))


command = "samtools mpileup -B -q 1 -Q 20  " + \
          " --positions " + snakemake.input.regions + \
          " -f " + snakemake.input.ref + \
          " " + snakemake.input.tumor + \
          " > " + snakemake.params.tumor_pileup + \
          " 2>> " + log_filename

f = open(log_filename, 'at')
f.write("## COMMAND: " + command + "\n")
f.close()
shell(command)

if snakemake.params.calling_type == "paired":
    command = "samtools mpileup -B -q 1 -Q 20  " + \
              " --positions " + snakemake.input.regions + \
              " -f " + snakemake.input.ref + \
              " " + snakemake.input.normal + \
              " > " + snakemake.params.normal_pileup + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)



# Output prefix


if snakemake.params.calling_type == "paired":
    prefix = os.path.splitext(os.path.splitext(snakemake.output.snp)[0])[0]
    command = "varscan somatic" \
              " " + snakemake.params.normal_pileup + \
              " " + snakemake.params.tumor_pileup + \
              " " + prefix + \
              " " + extra_params + \
              " --output-vcf" + \
              " --output-snp " + snakemake.output.snp + \
              " --output-indel " + snakemake.output.indel + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "rm " + snakemake.params.normal_pileup
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "rm " + snakemake.params.tumor_pileup
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

else:


    command = "varscan mpileup2snp " + snakemake.params.tumor_pileup + \
              " " + extra_params + \
              " --output-vcf" + \
              " > " + snakemake.params.snp + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "varscan mpileup2indel " + snakemake.params.tumor_pileup + \
              " " + extra_params + \
              " --output-vcf" + \
              " > " + snakemake.params.indel + \
              " 2>> " + log_filename

    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "rm " + snakemake.params.tumor_pileup
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    if os.stat(snakemake.params.snp).st_size > 0 or os.stat(snakemake.params.indel).st_size > 0:
        command = "Rscript " + os.path.abspath(os.path.dirname(__file__)) + "/combine_vcfs.R " + \
                  " " + snakemake.params.snp + \
                  " " + snakemake.params.indel + \
                  " " + snakemake.output.vcf + \
                  " >> " + log_filename + " 2>&1"

        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
        f.close()
        shell(command)

    else:
        command = "touch " + snakemake.output.vcf + " >> " + log_filename + " 2>&1"
        f = open(log_filename, 'at')
        f.write("## COMMAND: " + command + "\n")
        f.close()
        shell(command)

# # DELETE mileups
# command = "rm " + snakemake.params.tumor_pileup
# f = open(log_filename, 'at')
# f.write("## COMMAND: " + command + "\n")
# f.close()
# shell(command)

# command = "rm " + snakemake.params.normal_pileup
# f = open(log_filename, 'at')
# f.write("## COMMAND: " + command + "\n")
# f.close()
# shell(command)




# command = "varscan processSomatic" \
#           " " + prefix + ".snp.vcf" +\
#           " 2>> " + log_filename
#
# f = open(log_filename, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
#
# java -Xmx4G -jar /VarScan2.3.7.jar somaticFilter \
# /70fb9a87cec148389cb2f13615c12908/varcalling/VarScan2.snp.Somatic.hc.vcf \
# -indel-file /70fb9a87cec148389cb2f13615c12908/varcalling/VarScan2.indel.vcf \
# -output-file /70fb9a87cec148389cb2f13615c12908/varcalling/VarScan2.snp.Somatic.hc.filter.vcf