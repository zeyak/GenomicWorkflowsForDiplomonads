"  DATA
  PE = aa 236 240 run1_R1.fastq.gz run1_R2.fastq.gz
  PE = ab 229 232 run2_R1.fastq.gz run2_R2.fastq.gz
  PE = ac 235 238 run3_R1.fastq.gz run3_R2.fastq.gz
  PACBIO = pacbio_nano_concat.fastq.gz # Both of the long reads are given
  END

  PARAMETERS
  GRAPH_KMER_SIZE = auto
  USE_LINKING_MATES = 0
  CA_PARAMETERS =  cgwErrorRate=0.15
  KMER_COUNT_THRESHOLD = 1
  NUM_THREADS = 30
  JF_SIZE = 1380000000
  SOAP_ASSEMBLY = 0
  FLYE_ASSEMBLY = 1
  END"