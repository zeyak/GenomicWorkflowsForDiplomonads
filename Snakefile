configfile: "env/config.yaml"

rule all:
    input:
        expand("output/tmhmm/{pfam}.txt", pfam=config['pfam']),
        expand("output/signalp/{pfam}.gff",pfam=config['pfam'])

rule tmhmm:
    input:
        protein ="output/interproscan/processed_data/{pfam}.faa"
    output:
        "output/tmhmm/{pfam}.txt"
    conda:
         "env/ComparativeGenomicsAnalysis.yaml"
    script:
        "scripts/3_ComparativeGenomicsAnalysis/3_ProteinDomainLevel/tmhmm.py"

rule signalp:
    input:
        protein= "output/interproscan/processed_data/{pfam}.faa"
    output:
        "output/signalp/{pfam}.gff"
    conda:
         "env/ComparativeGenomicsAnalysis.yaml"
    script:
        "scripts/3_ComparativeGenomicsAnalysis/3_ProteinDomainLevel/signalp.py"

