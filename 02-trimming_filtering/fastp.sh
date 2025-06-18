#!/bin/bash
#Los archivos que acepta solo son archivos .fastq.gz 
#Author: Laura Barreales and Sara Lévano
#Start date: 14th May 2025
#Purpose: This script trims and filters raw sequences 

#Usage example: ./fastp.sh -m <MIN_QUAL> -l <MIN_LEN> -i <INPUT> -o <OUTPUT> -f <FASTQC_DIR> -c <MULTIQC_DIR>
#-v and -h are available to check script's version or ask for script's usage

#This directories should be the input so that this script runs correctly (in this context)
#input_dir="../00-raw_data/results"
#output_dir="./results"
#fastqc_dir="./results/fastqc"
#multiqc_dir="./results/multiqc"

#If user does not introduce any minimun quality or minimum length, minimum quality will be 20, and minimum length will be 40bp

#Link to fastp manual: https://open.bioqueue.org/home/knowledge/showKnowledge/sig/fastp 

#######################################################################

readonly version="Version: 1.3"

# Initializing empty logs
cat /dev/null > logs/*

#Usage
{ echo -e "Use -h for help or -v for the version. \nIntroudce input, output, fastqc and multiqc directories, but do not introduce more arguments unless you want to change filter parameters." 
 
    MIN_QUAL=20
    MIN_LEN=40 
    
    #Check if arguments are provided
    if [ $# -lt 4 ]; then
        echo -e "Need at least four arguments. \nRun $0 -h for usage help" >&2
        exit 1
    fi
    #Managing options and arguments
    while getopts "hvm:l:i:o:f:c:" opt; do
        case $opt in
            h) echo -e "Usage help:\nThis script filters and trims raw sequences using fastp."
            echo -e "Usage:\n\t$0 -m <MIN_QUAL> -l <MIN_LEN> -i <INPUT> -o <OUTPUT> -f <FASTQC_DIR> -c <MULTIQC_DIR>"
            exit 0;;
            v) echo "$version"
            exit 0;;
            m) MIN_QUAL="$OPTARG" ;;
            l) MIN_LEN="$OPTARG" ;;
            i) input_dir="$OPTARG" ;;
            o) output_dir="$OPTARG" ;;
            f) fastqc_dir="$OPTARG" ;;
            c) multiqc_dir="$OPTARG" ;;
            \?) echo -e "Invalid option\nUsage: $0 ..." >&2
                exit 1;; 
        esac
    done
    echo "Using minimun quality = $MIN_QUAL and minimun length = $MIN_LEN"


    # Checking if directory exists. If output_dir does not exist, it will be created
    echo -e "\nChecking if input directory exists" 
    if [ ! -d "$input_dir" ]
    then 
        echo "$input_dir does not path exist. Creating..."
        mkdir -p "$input_dir"
        echo "Created"
    else
        echo " "$input_dir" exists "
    fi
    
    echo -e "\nChecking if output directory exists" 
    if [ ! -d "$output_dir" ]
    then 
        echo "$output_dir does not path exist. Creating..."
        mkdir -p "$output_dir"
        echo "Created"
    else
        echo " "$out_dir" exists "
    fi
    
    echo -e "\nChecking if fastqc directory exists" 
    if [ ! -d "$fastqc_dir" ]
    then 
        echo "$fastqc_dir does not path exis. Creating..."
        mkdir -p "$fastqc_dir"
        echo "Created"
    else
        echo " "$fastqc_dir" exists "
    fi
    
    echo -e "\nChecking if multiqc directory exists" 
    if [ ! -d "$multiqc_dir" ]
    then 
        echo "$multiqc_dir does not path exist. Creating"
        mkdir -p "$multiqc_dir"
        echo "Created"
    else
        echo " "$multiqc_dir" exists "
    fi
} 2>> >(tee -a logs/all_files.err) >> >(tee -a logs/all_files.out)

# Processing fastqc files. As they are paired sequences they will be processed two by two
for file in "$input_dir"/*_1.fastq.gz; do 
    sample=$(basename "$file" "_1.fastq.gz") 
    file1="${input_dir}/${sample}_1.fastq.gz"
    file2="${input_dir}/${sample}_2.fastq.gz" 
    #Checking emptiness of the file
    { echo "Checking if $sample exists and is not empty"
    if [ -s "$file1" ] && [ -s "$file2" ]; then
        echo "Files exist and are not empty"
    else
        echo "One or both paired files are empty or missing" >&2
        exit 1
    fi
    #Checking if file has a .fastq.gz extension 
    echo "Checking whether file extension is .fastq.gz"
    if [[ "$file1" == *.fastq.gz && "$file2" == *.fastq.gz ]]; then 
        echo "File extensions are correct"
    else
        echo "Incorrect file extensions" >&2
        exit 1
    fi
    #Checking quality with fastqc
    echo "Checking read permissions"
    for f in "$file1" "$file2"; do
        if [[ ! -r "$f" ]]; then 
            echo "File $f not readable. Attempting to fix." >&2
            chmod +r "$f"
        fi
    done
    } 2>> >(tee -a logs/${sample}.err) >> >(tee -a logs/${sample}.out)
    #Initiating trimming
    echo "Processing $sample..."

    {
    fastp \
        -i "$file1" -I "$file2" \
        -o "$output_dir/${sample}_1_filtered.fastq.gz" \
        -O "$output_dir/${sample}_2_filtered.fastq.gz" \
        -q "$MIN_QUAL" -l "$MIN_LEN" \
        --detect_adapter_for_pe \
        --thread 4 \
        --html "$output_dir/${sample}_fastp.html" \
        --json "$output_dir/${sample}_fastp.json" \
        --trim_poly_x --poly_x_min_len 10 \
        --trim_front1 10 --trim_front2 10\
        --report_title "${sample}_fastp_report"  
        #i: input files from de input_dir
        #o: outputfiles that will be stored in the output directory
        #q: minimum quality 
        #l: minimum length after trimming
        #Detects adapters and trims them
        #Uses 4 threads
        #Generates .html and .json reports
        #Detects polyA and polyT that are minimum 10 bases long and trims them 
        #Each sample has its report 
    } 2>> >(tee -a logs/${sample}_fastp.err) >> >(tee -a logs/${sample}_fastp.out)
#Checking quality of filtered and trimmed sequences using fastqc and multiqc. REsults will be stored in fastqc_dir and multiqc_dir
    {
    echo "Executing FastQC for $sample" 
    fastqc "$output_dir/${sample}_1_filtered.fastq.gz" "$output_dir/${sample}_2_filtered.fastq.gz" -o "$fastqc_dir"
    } 2>> >(tee -a logs/fastqc/${sample}_fastqc.err) >> >(tee -a logs/fastqc/${sample}_fastqc.out)
done

{ echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$fastqc_dir" -n "multiqc_fastp_analysis.html" -o "$multiqc_dir"
} 2>> >(tee -a logs/multiqc/multiqc.err) >> >(tee -a logs/multiqc/multiqc.out)

echo -e "Analysis completed." >> >(tee -a logs/all_files.out)
