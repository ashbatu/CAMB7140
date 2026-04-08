#!/bin/bash
sra_accession="SRR8668771"
fastq_file="fastq/${sra_accession}.fastq.gz"
downsampled_fastq_file="downsampled_fastq/${sra_accession}.fastq.gz"

echo "Downsampling $sra_accession"
seqtk sample -s100 "$fastq_file" 1000000 | pigz > "$downsampled_fastq_file"
