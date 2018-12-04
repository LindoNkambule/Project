#!usr/bin/env bash
sudo mkdir /usr/local/pipeline/BWA_Alignment_Files

#Create a BWA index in the genomic reference
echo "==================================== INDEXING ==============================================="
bwa index GRCh38_latest_genomic.fna
echo "================================== DONE INDEXING ============================================"

for filename in *.fastq; do
  FastqFilename="$filename"
	sudo mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}

	samFilename=${FastqFilename//'.sam'/'.fastq'}
	sudo mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/SAM_Files

	bamFilename=${samFilename//'.bam'/'.sam'}
	sudo mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/BAM_Files

	sortedBAMfilename=${bamFilename//'.sorted.bam'/'.bam'}
	sudo mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/Sorted_BAM_Files

  statsFilename=${FastqFilename//'.stats'/.'fastq'}
	sudo mkdir /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/Stats

#Once the indexing is ready you can align the reads in the input file against the genomic reference:

#Single-end alignment (using BWA-MEM algorithm)
  bwa mem GRCh38_latest_genomic.fna $FastqFilename  > $samFilename (#this is the output)

#Converting SAM to BAM with samtools
  samtools view -S -b $samFilename > $bamFilename #We must specify our input -S for SAM format, then also specify the type of output we want -b  for BAM format

#Sorting our BAM file before we call variants
  samtools sort $bamFilename -o $sortedBAMfilename #The -o is for output

#Getting the basic stats for our BAM files
  bamtools stats -in $sortedBAMfilename > $statsFilename #The -in is for specifying the output

  mv *.sam /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/SAM_Files
  mv *.bam /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/BAM_Files
  mv *.sorted.bam /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/Sorted_BAM_Files
  mv *.stats /usr/local/pipeline/BWA_Alignment_Files/${filename%.*.*}/Stats
done

echo "=================================== ALIGNMENT DONE =========================================="
