""" Implement the following on python"""

"""
>awk '!seen[$1]++' prok_contamination_fmt6.txt | sort > prok_uniqe.txt
>awk '{print $1}' prok_uniqe.txt | perl /home/alejandro/scripts/extract_seq_2.0.pl assembly.fasta - > prok_uniqe.fasta
>perl /home/alejandro/scripts/extract_seq_not_equal.pl assembly.fasta filter1.txt > filter_clean.fasta #extract_seq_not_equal
>perl /home/alejandro/scripts/extract_seq_2.0.pl assembly.fasta filter1.txt > filter.fasta #extract_seq_2.0
>perl /home/lisa/Scripts/GC_count.pl filter.fasta
"""
from snakemake.shell import shell

putative_contamination= snakemake.input.putative_contamination
contamination = snakemake.output.contamination

shell(f"ContamintationFilter.py -i {putative_contamination} -o {contamination}")