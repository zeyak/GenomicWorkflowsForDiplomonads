configfile: "config.yaml"

rule all:
    input:
        (
            #expand("output/tmhmm/{pfam}.txt", pfam=config['pfam']),
            #expand("output/signalp/{pfam}.gff",pfam=config['pfam'])
            expand("output/spiro/fastqc/{sample}.html",sample=config['sample']),
            "output/spiro/multiqc/"
        )

rule tmhmm:
    input:
        protein="output/interproscan/processed_data/{pfam}.faa"
    output:
        "output/tmhmm/{pfam}.txt"
    conda:
        "envs/comparative_genomics.yaml"
    script:
        "scripts/comparative_genomics/protein_domain_level/run_tmhmm.py"

rule signalp:
    input:
        protein="output/interproscan/processed_data/{pfam}.faa"
    output:
        "output/signalp/{pfam}.gff"
    conda:
        "envs/comparative_genomics.yaml"
    script:
        "scripts/comparative_genomics/protein_domain_level/run_signalp.py"

rule fastqc:
    input:
        "resources/spiro/{sample}.fastq.gz"
    output:
        "output/spiro/fastqc/{sample}.html"
    conda:
        "envs/hybrid_genome_assembly.yaml"
    script:
        "scripts/hybrid_genome_assembly/read_preprocessing/check_read_quality.py"

rule multiqc:
    input:
        "resources/spiro/SRR8895272.fastq.gz"
    output:
        directory("output/spiro/multiqc/")
    conda:
        "envs/hybrid_genome_assembly.yaml"
    script:
        "scripts/hybrid_genome_assembly/read_preprocessing/get_multiqc.py"

rule makeblastdb:
    input:
        "resources/{type}/db/{db}.fasta"
    output:
        multiext("output/{type}/db/{db}",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto")
    params:
        outname="output/{type}/db/{db}"
    conda:
        "envs/blast.yaml"
    shell:
        'makeblastdb -dbtype nucl -in {input} -out {params.outname}'

rule blastn:
    input:
        query="resources/{type}/query/{query}.fasta",
        db="output/{type}/db/{db}.ndb"
    output:
        'output/{type}/{db}/{query}.blastn'
    params:
        perc_identity=95,
        outfmt=6,
        num_threads=30,
        max_target_seqs=1,
        max_hsps=1,
        db_prefix="output/{type}/db/{db}"
    conda:
        "envs/blast.yaml"
    script:
        "scripts/blastn.py"

rule setup_nr_db:
    output:
        protected(directory("/data/zeynep/databases"))
    conda:
        "envs/blast.yaml"
    script:
        "scripts/setup_nr_db.py"

rule blastdbcmd:
    input:
        "output/{type}/{db}/{prefix}"
    output:
        "output/{type}/{db}/{prefix}.blast.fasta"
    params:
        db_prefix="/data/zeynep/databases/nr"
    conda:
        "envs/blast.yaml"
    shell:
        "blastdbcmd -db {params.db_prefix} -entry_batch {input} > {output}"

rule scatter_fasta:
    input:
        "output/{type}/{db}/{prefix}.blast.fasta"
    output:
        "output/{type}/{db}/{prefix}_{n}.blast.fasta"
    params:
        n_partitions=10,
        n_partition="{n}"
    conda:
        "envs/blast.yaml"
    script:
        "scripts/scatter_fasta.py"
