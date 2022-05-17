# STAR mapping 
Working directory: 
```sh
/home/kukoll1/projects/rrg-ben/kukoll1/dmrt1/
```
# copying files with rsync
```sh
rsync -axvH --no-g --no-p ../../../for_lindsey/scripts/2022* .
```
# download the genome and annotation 
```sh
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
wget http://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_GCF_XBmodels.gtf
```
# uncompress the reference 
```sh
gunzip XENLA_10.1_genome.fa.gz
```

# index the genome 
```sh
#!/bin/sh
#SBATCH --job-name=star_index
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=5:00:00
#SBATCH --mem=16gb
#SBATCH --output=star_index.%J.out
#SBATCH --error=star_index.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2022_align_paired_fq_to_ref_for_star.sh path name_of_ref_fasta name_of_ref name_of_gtf_file

module load StdEnv/2020 star/2.7.9a

STAR --runThreadN 4 \
--runMode genomeGenerate \
--genomeDir ${1} \
--genomeFastaFiles ${1}/${2} \
--sjdbGTFfile ${1}/${3} \
--sjdbOverhang 99
```
# mapping data for each individual to the reference
For concatenated (dmrt1L) files: 
```sh
#!/bin/sh
#SBATCH --job-name=star_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=120:00:00
#SBATCH --mem=48gb
#SBATCH --output=star_align.%J.out
#SBATCH --error=star_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2022_align_paired_fq_to_ref_for_star.sh pathandname_of_ref path_to_paired_fq_filez prefix
# dmrt1L_55_R1_trim_cat.fastq.gz



module load StdEnv/2020 star/2.7.9a

STAR --genomeDir ${1} \
             --runThreadN 6 \
             --readFilesIn ${2}/${3}_R1_trim_cat.fastq.gz ${2}/${3}_R2_trim_cat.fastq.gz \
             --outFileNamePrefix ${3} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard

```
For non-concatenated (dmrt1L) and dmrt1S files:
``` sh
#!/bin/sh
#SBATCH --job-name=star_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=120:00:00
#SBATCH --mem=48gb
#SBATCH --output=star_align.%J.out
#SBATCH --error=star_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2022_align_paired_fq_to_ref_for_star.sh pathandname_of_ref path_to_paired_fq_filez prefix
# dmrt1L_55_R1_trim_cat.fastq.gz



module load StdEnv/2020 star/2.7.9a

STAR --genomeDir ${1} \
             --runThreadN 6 \
             --readFilesIn ${2}/${3}_R1_trim.fastq.gz ${2}/${3}_R2_trim.fastq.gz \
             --outFileNamePrefix ${3} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard
             
 ```
 Sort and Readgroups (Samtools)
 ```
 #!/bin/sh
#SBATCH --job-name=sort_and_rg
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=sort_and_rg.%J.out
#SBATCH --error=sort_and_rg.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2022_samtools_sort_and_readgroups.sh bamfile_prefix samplename

module load StdEnv/2020 samtools/1.12
# Sort both alignments
samtools sort ${1}.bam -o ${1}.sorted.bam

# add readgroups
module load picard/2.23.3

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${1}.sorted.bam O=${1}.sorted_rg.bam RGID=4 RGLB=${2} RGPL=ILLUMINA RGPU=${2} RGSM=${2}

# index
samtools index ${1}.sorted_rg.bam
```

