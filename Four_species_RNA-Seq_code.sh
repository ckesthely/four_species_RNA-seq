#!/bin/bash

# Name of the job
#SBATCH --job-name=Pipeline_4_bug_community

# Number of compute nodes
#SBATCH --nodes=2

# Number of cores, in this case one
#SBATCH --ntasks-per-node=8

# Walltime (job duration)
#SBATCH --time=60:00:00

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL

hostname
date
sleep 60

cd /dartfs-hpc/scratch/f006fcn/redo_4_bug/results/raw_data


#run fastqc on data
module load fastqc
fastqc *.fastq.gz

#make a new diRevtory and put the fastqc files into it
mkdir ../fastqc_results
mv *fastqc* ../fastqc_results
cd ../fastqc_results

eval "$(conda shell.bash hook)"
conda activate bioinfo-notebook
#run multiqc on the fastqc files
multiqc .

mkdir ../trim
cd ../trim

###read trimming and cutadapt
ls ../raw_data/*001.fastq.gz | while read x; do \
        sample=`echo "$x"`
        sample=`echo "$sample" | cut -d"/" -f3 | rev | cut -d"_" -f3- | rev`
        echo processing "$sample"
        
        cutadapt \
        -o ${sample}_R1.trim.001.fastq.gz \
        -p ${sample}_R2.trim.001.fastq.gz \
        ../raw_data/${sample}_R1_001.fastq.gz ../raw_data/${sample}_R2_001.fastq.gz \
        -m 10 -q 20 --nextseq-trim=10 > $sample.cutadapt.report
done


#arguments that can be added to the cutadapt protocol
###		-a #insert adapter here - 3' end adapter
###		-b #insert adapter here - either end adapter - don't use if you know which end the adapter is on
###		-g #insert adapter here - 5' end adapter
		

#print reports from cutadapt
ls *.cutadapt.report | while read x; do
	#print empty line
	echo `` \
	#print report name
	echo Printing $x
	#print empty line
	echo `` \
	#print report
	cat $x
done

for file in *.gz; do
	mv -- "$file" "${file/_S*_R/_R}"
	done


#read mapping
mkdir -p PA14/{Fwd,Rev}
mkdir -p Staph/{Fwd,Rev}
mkdir -p Strep/{Fwd,Rev}
mkdir -p Prev/{Fwd,Rev}
mkdir Reports

mv *.report Reports

cp Mix*R1*.gz PA14/Fwd
cp Mix*R2*.gz PA14/Rev
mv PA14*R1*.gz PA14/Fwd
mv PA14*R2*.gz PA14/Rev
mv lasR*R1*.gz PA14/Fwd
mv lasR*R2*.gz PA14/Rev

cp Mix*R1*.gz Staph/Fwd
cp Mix*R2*.gz Staph/Rev
mv Staph*R1*.gz Staph/Fwd
mv Staph*R2*.gz Staph/Rev

cp Mix*R1*.gz Strep/Fwd
cp Mix*R2*.gz Strep/Rev
mv Strep*R1*.gz Strep/Fwd
mv Strep*R2*.gz Strep/Rev

mv Mix*R1*.gz Prev/Fwd
mv Mix*R2*.gz Prev/Rev

mkdir -p ../alignment/{PA14,Staph,Strep,Prev}
cd ../alignment/PA14


#align reads to the genome using previously indexed genomes

ls ../../trim/PA14/Fwd/*.001.fastq.gz | while read x; do
  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../../trim/PA14/' and '_R')
  sample=`echo "$sample" | cut -d"/" -f6 | rev | cut -d"_" -f2- | rev`
  echo processing "$sample"

  # run bowtie2 for each sample
  bowtie2 -x /dartfs-hpc/rc/home/n/f006fcn/PA14/PA14 \
  	-1 ../../trim/PA14/Fwd/${sample}_R1.trim.001.fastq.gz -2 ../../trim/PA14/Rev/${sample}_R2.trim.001.fastq.gz \
  	-S ./${sample}.sam
done

cd ../Staph
ls ../../trim/Staph/Fwd/*.001.fastq.gz | while read x; do
  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../../trim/PA14/' and '_R')
  sample=`echo "$sample" | cut -d"/" -f6 | rev | cut -d"_" -f2- | rev`
  echo processing "$sample"

  # run bowtie2 for each sample
  bowtie2 -x /dartfs-hpc/rc/home/n/f006fcn/Staph/Staph \
  	-1 ../../trim/Staph/Fwd/${sample}_R1.trim.001.fastq.gz -2 ../../trim/Staph/Rev/${sample}_R2.trim.001.fastq.gz \
  	-S ./${sample}.sam
done

cd ../Strep
ls ../../trim/Strep/Fwd/*.001.fastq.gz | while read x; do
  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../../trim/PA14/' and '_R')
  sample=`echo "$sample" | cut -d"/" -f6 | rev | cut -d"_" -f2- | rev`
  echo processing "$sample"

  # run bowtie2 for each sample
  bowtie2 -x /dartfs-hpc/rc/home/n/f006fcn/Strep/Strep \
  	-1 ../../trim/Strep/Fwd/${sample}_R1.trim.001.fastq.gz -2 ../../trim/Strep/Rev/${sample}_R2.trim.001.fastq.gz \
  	-S ./${sample}.sam
done

cd ../Prev
ls ../../trim/Prev/Fwd/*.001.fastq.gz | while read x; do
  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../../trim/PA14/' and '_R')
  sample=`echo "$sample" | cut -d"/" -f6 | rev | cut -d"_" -f2- | rev`
  echo processing "$sample"

  # run bowtie2 for each sample
  bowtie2 -x /dartfs-hpc/rc/home/n/f006fcn/Prev/Prev \
  	-1 ../../trim/Prev/Fwd/${sample}_R1.trim.001.fastq.gz -2 ../../trim/Prev/Rev/${sample}_R2.trim.001.fastq.gz \
  	-S ./${sample}.sam
done



#sort and convert sam to bam
cd ../PA14

ls ./*.sam | while read x; do
    # save the file name
    sample=`echo "$x"`
    # get everything in file name before "/" (to remove './')
    sample=`echo "$sample" | cut -d"/" -f2 | cut -d"." -f1`
    echo processing "$sample"
                
    samtools sort ./${sample}.sam > ./${sample}_PA14.sorted.sam
done


cd ../Staph

ls ./*.sam | while read x; do
    # save the file name
    sample=`echo "$x"`
    # get everything in file name before "/" (to remove './')
    sample=`echo "$sample" | cut -d"/" -f2 | cut -d"." -f1`
    echo processing "$sample"
                
    samtools sort ./${sample}.sam > ./${sample}_Staph.sorted.sam
done


cd ../Strep

ls ./*.sam | while read x; do
    # save the file name
    sample=`echo "$x"`
    # get everything in file name before "/" (to remove './')
    sample=`echo "$sample" | cut -d"/" -f2 | cut -d"." -f1`
    echo processing "$sample"
                
    samtools sort ./${sample}.sam > ./${sample}_Strep.sorted.sam
done


cd ../Prev

ls ./*.sam | while read x; do
    # save the file name
    sample=`echo "$x"`
    # get everything in file name before "/" (to remove './')
    sample=`echo "$sample" | cut -d"/" -f2 | cut -d"." -f1`
    echo processing "$sample"
                
    samtools sort ./${sample}.sam > ./${sample}_Prev.sorted.sam
done


#make feature counts tables

#PA14
cd ../../
mkdir -p featureCounts/{PA14,Staph,Strep,Prev}
cd featureCounts/PA14

ls ../../alignment/PA14/*.sorted.sam | while read x; do
	# save the file name
	sample=`echo "$x"`
	# get everything in file name before "/" (to remove './')
	sample=`echo "$sample" | cut -d"/" -f5 | cut -d"." -f1`
	echo processing "$sample"
	
featureCounts \
	-a ../../gtf_files/PA14.gtf \
	-o ./${sample}_featureCounts_output.txt \
	../../alignment/PA14/${sample}.sorted.sam \
	-p -O -t CDS -g gene_id -f -d 10 --verbose
done

# set up an array to fill with shorthand sample names
myarray=()

# loop over htseq.counts files and extract 2nd column (the raw read counts) using 'cut' command
while read x;  do
  # split up sample names to remove everything after "-"
  sname=`echo "$x"`
  sname=`echo "$sname" | rev | cut -d"_" -f3- | rev`
  # extract second column of file to get read counts only
  echo counts for "$sname" being extracted
  cut -f7 $x > "$sname".tmp.counts
  # save shorthand sample names into an array  
  sname2="$sname"
  myarray+=($sname2)
done < <(ls -1 *output.txt | sort)

#extract gene IDs and gene names
cut -f1 Mix_1_PA14_featureCounts_output.txt > genes.txt

# paste gene IDs into a file with each to make the gene expression matrix
paste genes.txt *.tmp.counts > tmp_all_counts.txt

# look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > col_names.txt
cat col_names.txt

# add 'gene_name' to colnames
cat <(echo "ENSEMBL_ID") <(echo "gene_name") col_names.txt > col_names_full.txt
cat col_names_full.txt

# make a file to fill
touch all_counts.txt

# use the 'cat' command (concatenate) to put all tmp.counts.txt files into all_counts.txt
cat <(cat col_names_full.txt | sort | paste -s) tmp_all_counts.txt > PA14_all_counts.txt


#remove temp files
rm -f *tmp*

 
#Staph

cd ../Staph

ls ../../alignment/Staph/*.sorted.sam | while read x; do
	# save the file name
	sample=`echo "$x"`
	# get everything in file name before "/" (to remove './')
	sample=`echo "$sample" | cut -d"/" -f5 | cut -d"." -f1`
	echo processing "$sample"

featureCounts \
	-a ../../gtf_files/Staph.gff \
	-o ./${sample}_featureCounts_output.txt \
	../../alignment/Staph/${sample}.sorted.sam \
	-p -O -t CDS -g gene_id -f -d 10 --verbose
done

# set up an array to fill with shorthand sample names
myarray=()

# loop over htseq.counts files and extract 2nd column (the raw read counts) using 'cut' command
while read x;  do
  # split up sample names to remove everything after "-"
  sname=`echo "$x"`
  sname=`echo "$sname" | rev | cut -d"_" -f3- | rev`
  # extract second column of file to get read counts only
  echo counts for "$sname" being extracted
  cut -f7 $x > "$sname".tmp.counts
  # save shorthand sample names into an array  
  sname2="$sname"
  myarray+=($sname2)
done < <(ls -1 *output.txt | sort)

#extract gene IDs and gene names
cut -f1 Mix_1_Staph_featureCounts_output.txt > genes.txt

# paste gene IDs into a file with each to make the gene expression matrix
paste genes.txt *.tmp.counts > tmp_all_counts.txt

# look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > col_names.txt
cat col_names.txt

# add 'gene_name' to colnames
cat <(echo "ENSEMBL_ID") <(echo "gene_name") col_names.txt > col_names_full.txt
cat col_names_full.txt

# make a file to fill
touch all_counts.txt

# use the 'cat' command (concatenate) to put all tmp.counts.txt files into all_counts.txt
cat <(cat col_names_full.txt | sort | paste -s) tmp_all_counts.txt > Staph_all_counts.txt


#remove temp files
rm -f *tmp*


#Strep

cd ../Strep

ls ../../alignment/Strep/*.sorted.sam | while read x; do
	# save the file name
	sample=`echo "$x"`
	# get everything in file name before "/" (to remove './')
	sample=`echo "$sample" | cut -d"/" -f5 | cut -d"." -f1`
	echo processing "$sample"

featureCounts \
	-a ../../gtf_files/Strep.gtf \
	-o ./${sample}_featureCounts_output.txt \
	../../alignment/Strep/${sample}.sorted.sam \
	-p -O -t CDS -g gene_id -f -d 10 --verbose
done

# set up an array to fill with shorthand sample names
myarray=()

# loop over htseq.counts files and extract 2nd column (the raw read counts) using 'cut' command
while read x;  do
  # split up sample names to remove everything after "-"
  sname=`echo "$x"`
  sname=`echo "$sname" | rev | cut -d"_" -f3- | rev`
  # extract second column of file to get read counts only
  echo counts for "$sname" being extracted
  cut -f7 $x > "$sname".tmp.counts
  # save shorthand sample names into an array  
  sname2="$sname"
  myarray+=($sname2)
done < <(ls -1 *output.txt | sort)

#extract gene IDs and gene names
cut -f1 Mix_1_Strep_featureCounts_output.txt > genes.txt

# paste gene IDs into a file with each to make the gene expression matrix
paste genes.txt *.tmp.counts > tmp_all_counts.txt

# look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > col_names.txt
cat col_names.txt

# add 'gene_name' to colnames
cat <(echo "ENSEMBL_ID") <(echo "gene_name") col_names.txt > col_names_full.txt
cat col_names_full.txt

# make a file to fill
touch all_counts.txt

# use the 'cat' command (concatenate) to put all tmp.counts.txt files into all_counts.txt
cat <(cat col_names_full.txt | sort | paste -s) tmp_all_counts.txt > Strep_all_counts.txt


#remove temp files
rm -f *tmp*


#Prev

cd ../Prev

ls ../../alignment/Prev/*.sorted.sam | while read x; do
	# save the file name
	sample=`echo "$x"`
	# get everything in file name before "/" (to remove './')
	sample=`echo "$sample" | cut -d"/" -f5 | cut -d"." -f1`
	echo processing "$sample"

featureCounts \
	-a ../gtf_files/Prev.gtf \
	-o ./${sample}_featureCounts_output.txt \
	../../alignment/Prev/${sample}.sorted.sam \
	-p -O -t CDS -g gene_id -f -d 10 --verbose
done

# set up an array to fill with shorthand sample names
myarray=()

# loop over htseq.counts files and extract 2nd column (the raw read counts) using 'cut' command
while read x;  do
  # split up sample names to remove everything after "-"
  sname=`echo "$x"`
  sname=`echo "$sname" | rev | cut -d"_" -f3- | rev`
  # extract second column of file to get read counts only
  echo counts for "$sname" being extracted
  cut -f7 $x > "$sname".tmp.counts
  # save shorthand sample names into an array  
  sname2="$sname"
  myarray+=($sname2)
done < <(ls -1 *output.txt | sort)

#extract gene IDs and gene names
cut -f1 Mix_1_Prev_featureCounts_output.txt > genes.txt

# paste gene IDs into a file with each to make the gene expression matrix
paste genes.txt *.tmp.counts > tmp_all_counts.txt

# look at the contents of the array we made with shorthand sample names
echo ${myarray[@]}

# print contents of array into text file with each element on a new line
printf "%s\n" "${myarray[@]}" > col_names.txt
cat col_names.txt

# add 'gene_name' to colnames
cat <(echo "ENSEMBL_ID") <(echo "gene_name") col_names.txt > col_names_full.txt
cat col_names_full.txt

# make a file to fill
touch all_counts.txt

# use the 'cat' command (concatenate) to put all tmp.counts.txt files into all_counts.txt
cat <(cat col_names_full.txt | sort | paste -s) tmp_all_counts.txt > Prev_all_counts.txt


#remove temp files
rm -f *tmp*
