######################################
# wrapper for rule: MuSE
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## WRAPPER: MuSE \n##\n")
f.close()

shell.executable("/bin/bash")

# version = str(subprocess.Popen("vardict-java 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
# f = open(log_filename, 'a+')
# f.write("## VERSION: vardict-java "+version+"\n")
# f.close()

command = "awk '{{ print $1, $2, $3 }}' " + snakemake.input.regions + " > " + snakemake.params.threecol_bed_tmp

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "MuSE call " +\
            " -O "+ snakemake.params.dir +\
            " -f "+ snakemake.input.ref +\
            " -l "+ snakemake.params.threecol_bed_tmp +\
            " " + snakemake.input.tumor +\
            " "+ snakemake.input.normal +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm " + snakemake.params.threecol_bed_tmp

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "MuSE sump " +\
            " -I "+ snakemake.params.txt_res +\
            " -D "+ snakemake.input.dbsnp +\
            " -O "+ snakemake.output.vcf +\
            " -E " +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
