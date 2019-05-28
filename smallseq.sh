# -----Clean & Concatenate FASTQ files-------------

# Quality check
fastqc *.fq
multiqc .

# Identify adapters | bbmerge.sh is part of the BBtools package | dnapi.py is part of DNApi
bbmerge.sh in=Reads_1.fq in2=Reads_2.fq mininsert=17 outa=adapterst.fa
dnapi.py --show-all Reads.fq

# Trim the adapters | Paired-end | Single-end
cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACTACAGCATCTCGTATGCCGTCTTCTGCTTG -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -o Reads_1.trim.fq -p Reads_2.trim.fq Reads_1.fq Reads_2.fq
cutadapt -a TGGAATTCTCGGGTGCCAAGG -o Read_1.trim.fq Read_1.fq
umi_tools extract --stdin=Read_1.fq --log=output.log --stdout=Read_1.trim.fq --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'

# Paired-end | Generate reverse complement of Reads_2.trim.fq and concatenate 
seqkit seq Reads_2.trim.fq -r -p -o Reads_2.trimrc.fq
cat Reads_1.trim.fq Reads_2.trimrc.fq > merged.fq

# Create config txt file as follows | list.txt
#path to file	                three letter code
Merged1.fastq                   S01
Merged2.fastq                   S02
Merged3.fastq                   S03
Merged4.fastq                   S04


# ----Extracting human mirBase mature and precursor sequences-----------

# Download
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

# Downloaded files have white spaces that need to be removed (for miRDeep2)
zcat hairpin.fa.gz | cut -f1 -d" ">hairpin_noSpace.fa
zcat mature.fa.gz | cut -f1 -d" " >mature_noSpace.fa

# Remove non-canonical nucleotides and convert to SL | fastaparse.pl is part of the miRDeep 2 package | fasta_formatter is part of FASTX-Toolkit
fastaparse.pl hairpin_noSpace.fa -b | fasta_formatter -w 0  >hairpin_cleaned.fa
fastaparse.pl mature_noSpace.fa -b | fasta_formatter -w 0 >mature_cleaned.fa

# grep out the human miRNAscd ..
grep -A1 "hsa" mature_cleaned.fa  | grep -v -- "^--$" >hsa_matureMirnas.fa
grep -A1 "hsa" hairpin_cleaned.fa | grep -v -- "^--$" >hsa_hairpinMirnas.fa


# ----Indexing human genome-----------

# Remove whitespace and indexing | remove_white_space_in_id.pl is part of the miRDeep 2 package | bowtie-build is part of the bowtie package
remove_white_space_in_id.pl genome.fa > genome_rs.fa
bowtie-build genome_rs.fa /bowtieindex/removedspaces


# ----miRDeep2-------------------

#Begin miRDeep process with mapper.pl (with config file)
mapper.pl list.txt -d -e -m -l 18 -v -q -h -j -p /bowtieindex/removedspaces -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf

# Output | FASTA file with processed reads | ARF file with with mapped reads

# Identify conserved miRNAs and predict novel miRNAs with miRDeep2
miRDeep2.pl /genome/genome_rs.fa reads_collapsed_vs_genome.arf hsa_matureMirnas.fa none hsa_hairpinMirnas.fa -t hsa -r results -P -v -c 2>report.log

# Output | A spreadsheet and an HTML file with an overview of all detected miRNAs in the deep sequencing input data.
