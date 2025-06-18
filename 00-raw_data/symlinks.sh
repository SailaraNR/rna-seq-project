#!/bin/bash
#Author: Laura Barreales y Sara Lévano
#Start date: 9th may 2025
#Purpose: This script will create symbolic links in the current directory which will enable the user to access fasta files of our classmates
#Usage example: symlinks.sh -d absolute/path/to/directory

#Arguments expected: -d absolute/path/to/directory

############################################################################################

readonly version="version 2.0"

# Initializing empty logs
cat /dev/null > logs/*


# Explain the code
{ echo -e "Use: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" 

#Check if arguments are provided
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help" >&2
        exit 1
fi

#Managing options and arguments
while getopts ":d:hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script creates on symbolic link per file.fastq in the directory provided" 
                   echo -e "Usage: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" 
                   exit 0;;
                v) echo $version 
                   exit 0;;
                d) directory="$OPTARG"
                   echo "Proccessing directory: $directory ..." 
                    ;;
                :) echo "Option -$OPTARG requieres an argument" >&2 #if only the flag is written, the script exits
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" >&2
                   exit 1 ;;

        esac
done
} 2>> >(tee -a logs/all_samples.err) >> >(tee -a logs/all_samples.out) 

#Accesing fastafiles and creating symlinks named after the sample without file extension
for file in "$directory"/*; do
    echo -e "\n" } >> >(tee -a logs/${name}.out) 
    sample=$(basename "$file") 
    name=${sample%.*} 
    { echo "This is $sample" 
    #Checking emptiness of the file
    echo "Checking if the file exists and it's not empty" 
    if [ -s "$file" ]; then
            echo "File exist and is not empty"  
    else
            echo "file is empty" >&2
            exit 1
    fi 

    #Checking if file has a .fastq.gz extension. If its not, that file will be skkiped
    echo "Checking whether file extension is .fastq.gz" 
    if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" 
    else
            echo "File extension not correct. Must be .fastq.gz file. Skipping $sample " >&2
            continue 
    fi 

    #Cheking if the file readable and executable
    echo "Checking $sample permissions" 

    if [[ -r "$file" ]]; then
            echo "File is readable. We can continue" 
    else
            echo "File is not readable . Fixing" >&2
            chmod +rx "$file"
            echo "Fixed" >&2
    fi

    
    #Creating the links in the current directory
    echo "Creating symbolic link for $sample ..." 
    ln -s $file $sample   
    mv $sample results/
    #Checking if symlink creation was successful
    if [ $? != 0 ]; then
        echo "Something went wrong during symbolic links creation for file $sample" >&2
        exit 1  
    fi
    echo "Symbolic link created for $file and named $sample" 
    } 2>> >(tee -a logs/${name}.err) >> >(tee -a logs/${name}.out) 
done


echo -e "\n" >> >(tee -a logs/all_files.out)
echo "Symbolic links successfully created in $(pwd)" >> >(tee -a logs/all_samples.out)

