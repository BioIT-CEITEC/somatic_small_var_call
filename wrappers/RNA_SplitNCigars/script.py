######################################
# wrapper for rule: RNA_SplitNCigars
######################################
import os.path
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## wrapper: RNA_SplitNCigars \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("gatk SplitNCigarReads --version true 2>&1 | grep \"[Vv]ersion:\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: gatk "+version+"\n")
f.close()

command = "gatk --java-options \"-Xmx2g\" SplitNCigarReads" + \
               " --tmp-dir " + os.path.dirname(snakemake.output.bam) +\
               " -R "+ snakemake.input.ref +\
               " -I "+ snakemake.input.bam +\
               " -O " + snakemake.output.bam +\
               " >> " + log_filename + " 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = " mv " + snakemake.params.bai + " " + snakemake.output.bai

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# shell("rm " + snakemake.input.bam)
