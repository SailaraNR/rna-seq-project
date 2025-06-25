#!/bin/bash
#Author: Laura Barreales y Sara Lévano
#Start date: 1st may 2025
#Purpose: This scripts accepts a file.txt containing a list of accessions of the samples you want to download
#Usage example: 00-data.sh -f SRR_Acc_List.txt

#Arguments expected: -f file.txt
#Raw sequences will be downloaded in the current directory

############################################################################################
readonly version="Version 1.0"
# Initializing empty logs
cat /dev/null > logs/*.out
cat /dev/null > logs/*.err

# explain the code
{ echo "Use -f for the file, -h for help and -v for the version" 

#Check if arguments are provided
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help" >&2
        exit 1
fi

#Managing options and arguments
while getopts ":f:hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script downloads rawdata if a file.txt with SRR accessions is given" 
                   echo -e "Usage example:\n\t ./$0 -h: Provides help \n\t ./$0 -v: Tells script's version\n\t ./$0 -f file.txt: Downloads rawdata from de accessions in file.txt" 
                   exit 0;;
                v) echo $version 
                   exit 0;;
                f) echo "Abriendo el archivo -f $file ..."  
                   file="$OPTARG" ;;
                :) echo "Option -$OPTARG requieres an argument" >&2 #if only the flag is written, the script exits
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" >&2
                   exit 1 ;;

        esac
done
} 2>> >(tee -a logs/all_samples.err) >> >(tee -a logs/all_samples.out)

#Checking emptiness of the file
{ echo "Checking if the file exists and it's not empty" 
if [ -s $file ]; then
        echo "File exist and is not empty"  
else
        echo "file is empty"  >&2
        exit 1
fi 

#Checking if file has a .txt extension
echo "Checking whether file extension is .txt" 
if [ "$file" == *.txt ]; then
        echo "File extension is correct" 
else
        echo "File extension not correct. Must be .txt file" >&2
        exit 1
fi 

#Cheking if the file readable and executable
echo "Checking $file permissions" 

if [[ -r "$file" && -x "$file" ]]; then
        echo "File is readable and executable" 
elif [ -r "$file" ]; then
        echo "File is readable but not executable. Fixing." >&2
        chmod +x "$file"
        echo "fixed" >&2
elif [ -x "$file" ]; then
        echo "File is executable but not readable" >&2
        chmod +r "$file"
        echo "Fiexed" >&2
else
        echo "File is neither readable nor executable. Fixing" >&2
        chmod +rx "$file"
        echo "Fixed" >&2
fi 
} 2>> >(tee -a logs/all_samples.err) >> >(tee -a logs/all_samples.out)

#Accesions will be read one by one and the corresponding raw_seqs will be downloaded
#Converts .sra files into .fastq and it splits into two files if samples are paired
#.sra files will be removed as they won't be used anymore
{ echo "Processing accesions from $file..." 
while read -r ACCESSION; do
  name=$(sample%.*)
  [[ -z "$ACCESSION" ]] && continue 
  echo "Downloading $ACCESSION" 
  prefetch "$ACCESSION" 
  fasterq-dump --split-files "$ACCESSION"
  echo "Compressing ${ACCESSION}.fastq..."
  gzip "${ACCESSION}"*.fastq
done < "$file"

echo "Remoning .sra files..."
rm -r *.sra

echo "Finished. FASTQ files downloaded in $(pwd)" 
} 2>> >(tee -a logs/all_${name}.err) >> >(tee -a logs/all_${name}.out)

echo "Program finished."
