# !/bin/bash
#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
#Purpose: This scripts counts the number of read of an alignment

# Usage example: ./counts_alignment.sh -A <genome_annotation.gtf> -d <dir_sorted_bam>
#-v and -h are available for help and version

#This directories should be the input so that this script runs correctly (in this context)
#dir_sorted_bam="../04-raw_data/results"

#Output: There will be created a summary and a .txt file for each sample

# Manual to FeatureCounts: https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html
#######################################################################

readonly version="versión 1.3"

# Initializing empty logs
cat /dev/null > logs/

#Checking if .out and .err files exists, if not they will be created
while getopts "hvA:d:" opt; do
    case $opt in
        v) echo "$version" 
        exit 0;;
        h) echo -e "This script counts the number of reads of an aligment using ./$0 -A genome_annotation.gtf -d dir_sorted_bam"
        exit 0;;
        A) GTF="$OPTARG";;
        d) alignment="$OPTARG";;
        \?) echo -e "Invalid option: This script counts the number of reads of \n./$0 -A genome_annotation.gtf -d dir_sorted_bam"
        exit 1;;
    esac
done 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

# File validation: Check if the variables exist
# El GTF no es nuestro, así que tampoco podemos cambiar los permisos, solo verificar si existe o no
if [[ ! -e "$GTF" ]] || [[ ! -r "$GTF" ]]; then
        echo "Error: $GTF does not exist or does not have read permissions" >&2
        exit 1
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

#Comprobar si existe 04/results (en principio debería existir porque se lo damos predetermiando
#pero podemos mirar si está vacío)
{
    if [[ ! -e "$alignment" ]]; then
        echo "Error: $alignment does not exist" >&2
        exit 1
    fi

    if [[ ! -r "$alignment" ]]; then
        echo "$alignment is not readable. Fixing"
        chmod +r "$alignment"
        echo "Fixed"
    fi

    if [[ ! -x "$alignment" ]]; then
        echo "$alignment is not executable. Fixing"
        chmod +x "$alignment"
        echo "Fixed"
    fi
} 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

for file in ${alignment}/SRR*/*.bam; do # 04-reads_alignment/results en nuestro caso sería el $alignment
    sample=$(basename "$file" "_sorted.bam")
    {
    echo -e "\nChecking if $sample exists and it's not empty"
    if [[ -f "$file" && -s "$file" ]]; then
        echo "$sample exists and is not empty"
    else
        echo "$sample is not a file or is empty"  >&2
        exit 1
    fi

    echo "Checking $sample permissions"
    if [[ ! -r "$file" ]]; then
        echo "File does not have read. Fixing..." >&2
        chmod +r "$file" && echo "Permissions fixed" || echo "Could not change permisions" >&2
        continue
    fi

    #Create output file
    echo "Creating output $sample.txt file..."
    sample_out="results/$sample.txt"
    if ! [[ -e "$sample_out" ]]; then
            echo "Output file does not exists, creating..."
            touch $sample_out
    fi

    #Counting reads
    featureCounts -p -O -T 6 -a "$GTF" -o "$sample_out" "$file" && echo "Counts for "$sample" done" || echo "Counts for "$sample" failed"
    }  2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out)
done 

echo "Counting has finished"
# -p species that fragments (or templates) will be counted instead of reads. This is only applicable for pairs
# -t type of fragments (exon, intron...)
# -O assigns reads to all their overlapping meta-features.
# -T specifies the number (6) of threads to be used.
# -a is the genome annotation file ($GTF).
# -o specifies the name of the output file, which includes the read counts ($output).
# $file (loop) in $alignment is an alignment file: in this file, the reads we want to count are aligned with the GTF

