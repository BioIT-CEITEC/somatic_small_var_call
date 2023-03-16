######################################
# wrapper for rule: merge_variant_callers
######################################
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: merge_variant_callers \n##\n")
f.close()

def number_of_variants(fname):
    if os.stat(fname).st_size > 0:
        with open(fname) as f:
            j = 0
            for i, l in enumerate(f):
                if "#CHROM" in l:
                    j = i
        return i - j
    else:
        return 0

shell.executable("/bin/bash")

version = str(subprocess.Popen("gatk MergeVcfs --version true 2>&1 | grep \"[Vv]ersion:\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: gatk "+version+"\n")
f.close()

callers = list()
merge_variant_callers_text = ""
for input_vcf in snakemake.input.vcfs:
    if number_of_variants(input_vcf) > 0:
        callers.append(input_vcf)
        merge_variant_callers_text = merge_variant_callers_text + "-I " + input_vcf + " "

if len(callers) == 0:
    command = "cp "+snakemake.input.vcfs[0]+" "+snakemake.output.not_filtered_vcf
elif len(callers) == 1:
    command = "cp "+callers[0]+" "+snakemake.output.not_filtered_vcf
else:
    command = "gatk MergeVcfs " + \
               merge_variant_callers_text +\
               " -R "+ snakemake.input.ref +\
               " -D " + snakemake.input.dict +\
               " -O " + snakemake.output.not_filtered_vcf +\
               " >> " + log_filename + " 2>&1 "


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

# after_merge_processing
command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_after_merge.R "+\
            snakemake.output.not_filtered_vcf + " " +\
            snakemake.output.tsv + " " +\
            str(snakemake.params.min_variant_frequency) + " " +\
            str(snakemake.params.min_callers_threshold) + " " +\
            str(snakemake.params.min_var_reads_threshold) + " " +\
            snakemake.params.tmp_dir +\
            " >> " + log_filename + " 2>&1"

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)
