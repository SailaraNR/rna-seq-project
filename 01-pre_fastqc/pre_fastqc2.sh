#!/bin/bash
#Author: Laura Barreales and Sara Lévano
#Start date: 4th May 2025
#Purpose. This script use raw sequences to analyze its quality. Only accepts .fastq.gz files

#Usage example: pre-fastqc.sh -i <input_dir> -0 <output_dir>
#-v and -h are available to check script's version or ask for script's usage

#These must be the input paht and output path in order to correctely running the script
#input_path="../00-raw_data/results"
#output_path="./results"

#Output: For each sample there will be created two files: a .zip file and a .html file which has the analysis for that sample
#Another file called multiqc_analysis.html will be created. This file collects all de sample's analysis in one .html file
#All this files will be stored in a folder called resutls in the specified output directory
#Additionally each sample has two log files .out and .err where the stdou and stderr will be written
#There are two logfiles .out and .err, which refers not to each sample but to the first checking of characteristics

#Manual de FastQC: https://github.com/s-andrews/FastQC/blob/master/fastqc
#Manual de MultiQC. No es de github, pero he accedido desde github:https://docs.seqera.io/multiqc/getting_started/running_multiqc 

#######################################################################

readonly version="Version 2.0 "
# Initializing empty logs
cat /dev/null > logs/*


#explain the code
{ echo -e "Use -h for help or -v for the version. \nUse -i to introduce input directory and -o to introduce output directory for fastqc outputfiles" 


#Check if arguments are provided
if [ $# -eq 0 ]; then
        echo -e "Need at least one argument. \nRun $0 -h for usage help" >&2
        exit 1
fi

#Managing options and arguments
while getopts ":hvo:i:" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script checks the quality of the previously downloaded raw sequences" 
            echo -e "Usage example:\n\t./$0 -h: Provides help \n\t./$0 -v: Tells script's version\n\t./$0: checks quality score of raw sequences" \
            "\n\t./$0 -i /path/to/input/direcyoty \n\t./$0 -o /path/to/output/direcyoty" 
            exit 0;;
        v) echo "$version" 
            exit 0;;
        i) input_path="$OPTARG"
            echo "Procesando el directorio: $input_path ..." ;;
        o) output_path="$OPTARG"
            echo "Procesando el directorio: $output_path ..." ;;
        :) echo "Option -$OPTARG requieres an argument" >&2 #si pone solo -d no sirve
            exit 1 ;;    
        \?) echo -e "Invalid option\n\ 
            Usage example: $0 " >&2
            exit 1;; #si pone simbolos raros no sirve
    esac
done


# Checking if directory exists. If output_dir does not exist, it will be created
echo -e "\nChecking if output directory exists" 
if [ ! -d "$input_path" ]
then 
     echo "Input path exists"
else 
     if [ "$(ls $input_path)" ]
     then  
         echo "Directory is not empty"
     else 
         echo "Directory is empty" >&2
         exit 1 
     fi
fi

echo -e "\nChecking if output directory exists" 
if [ ! -d "$output_path" ]
    echo "Creating output directory as it does not exist"
    mkdir -p "$output_path"
then 
     echo "Input path exists"
else 
     if [ "$(ls "$output_path")" ]
     then  
         echo "Directory is not empty"
     else 
         echo "Directory is empty" >&2
     fi
fi
} 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out) 



for file in "$input_path"/*; do
        #Checking emptiness of the file
        sample=$(basename "$file") 
        name=${sample%.*} #output files name (does not have the extension)
        { echo -e "\nChecking if $sample exists and it's not empty" 
        if [[ -f "$file" && -s "$file" ]]; then
            echo "$sample exist and is not empty"  
         else
            echo "$sample is not a file or is empty"  >&2
            exit 1
        fi 

        #Checking if file has a .fastq.gz extension if it is not, the script will skip that file
        echo "Checking whether file extension is .fastq.gz" 
        if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" 
            #Cheking if the file readable
            echo "Checking $sample permissions" 
            if [[ ! -r "$file" ]]; then
                echo "File is not readable. Fixing..." >&2
                chmod +r "$file" && echo "Permissions fixed" || echo "Could not change permisions" >&2
                    continue
            fi
            
            #Checking quality with fastqc
            echo "Checking quality score in $sample" 
            fastqc $file -o $output_path   
            #Checking if the analysis was successful
            if  [ $? -eq 0 ]; then
                echo "Fastqc analysis completed for $sample" 
            else
                echo "Fastqc analysis failed for $sample" >&2
            fi    
        else
            echo "this is not .fastq.gz, moving to the next..." >&2
            continue

        fi;
        } 2>> >(tee -a logs/${name}.err) >> >(tee -a logs/${name}.out) 
done 

{ echo -e "\nNow Multiqc will run in order to make a summary of the analysis" 
    echo "But first it will check if files' extensions are .html after fastqc analysis" 

    #Using multiqc to collect all the analysis in one .html file called: multiqc_analysis.html
    multiqc $output_path -n multiqc_analysis.html 
    echo -e "Multiqc has finished its job. \nYou can now check for an .html file in $(pwd) \nto find the analysis" 
} 2>> >(tee -a logs/multiqc.err)  >> >(tee -a logs/multiqc.out) 


