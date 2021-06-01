######################################
# wrapper for rule: merge_variants_in_samples
######################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_variants_in_samples \n##\n")
f.close()

shell.executable("/bin/bash")


command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/merge_variants_in_samples.R "+\
            snakemake.output.tsv_for_vep + " " +\
            " ".join(snakemake.input.var_tabs) + " " +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)





# # after_merge_processing
# command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_after_merge.R " + not_filtered_vcf_file + " " + snakemake.output.tsv + " " + snakemake.wildcards.subclass
#
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
#
# shell(command)
