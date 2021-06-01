######################################
# wrapper for rule: postprocess_somaticseq_variants
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## WRAPPER: postprocess_somaticseq_variants \n##\n")
f.close()

shell.executable("/bin/bash")

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/postprocess_somaticseq_variants.R " +\
            snakemake.input.snv + " " +\
            snakemake.input.indel + " " +\
            snakemake.output.var_tab + " " +\
            snakemake.params.calling_type + " " +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
