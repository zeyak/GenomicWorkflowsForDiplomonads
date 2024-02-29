configfile: "env/config.yaml"

rule all:
    input:
        #expand("output/tmhmm/{pfam}.txt", pfam=config['pfam']),
        #expand("output/signalp/{pfam}.gff",pfam=config['pfam'])
        expand("output/spiro/fastqc/{sample}.html", sample=config['sample']),
        "output/spiro/multiqc/"

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

rule fastqc:
    input:
        "resources/spiro/{sample}.fastq.gz"
    output:
        "output/spiro/fastqc/{sample}.html"
    conda:
         "env/HybridGenomeAssemblyWorkflow.yaml"
    script:
        "scripts/1_HybridGenomeAssemblyWorkflow/1_ReadsPreprocessing/ReadQualityCheck.py"

rule multiqc:
    input:
        "resources/spiro/SRR8895272.fastq.gz"
    output:
        directory("output/spiro/multiqc/")
    conda:
         "env/HybridGenomeAssemblyWorkflow.yaml"
    script:
        "scripts/1_HybridGenomeAssemblyWorkflow/1_ReadsPreprocessing/MultiqcReport.py"
