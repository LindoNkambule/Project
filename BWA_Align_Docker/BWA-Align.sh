#!usr/bin/env bash
mkdir /usr/local/pipeline/BWA_Alignment_Files
#cd /usr/local/pipeline/BWA_Alignment_Files
#wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz    #Download the latest humnan reference genome
#gunzip -k GRCh38_latest_genomic.fna.gz                                # Unzip the file gz file using gunzip and keep (-k) both the decompressed file and orginal gz (compressed) file

#Create a BWA index in the genomic reference
#echo "==================================== INDEXING ==============================================="
#bwa index -a bwtsw GRCh38_latest_genomic.fna
#echo "================================== DONE INDEXING ============================================"

for filename in *.fastq; do
  FastqFilename="$filename"
	mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*}

#Once the indexing is ready you can align the reads in the input file against the genomic reference:

#Single-end alignment (using BWA-MEM algorithm)
  bwa mem GRCh38_latest_genomic.fna $FastqFilename  > output.sam #this is the output
#Converting SAM to BAM with samtools
  samtools view -S -b output.sam > output.bam #We must specify our input -S for SAM format, then also specify the type of output we want -b  for BAM format

#Sorting our BAM file before we call variants
  samtools sort output.bam -o output.sorted.bam #The -o is for output

  sudo mv *.sam *.bam /usr/local/pipeline/BWA_Alignment_Files/${filename%.*}
done

echo "=================================== ALIGNMENT DONE =========================================="
