# !/bin/bash
#En este script se crearan los links simbólicos para poder acceder a los archivos fasta de nuestas
#compañeras de clase.
#Author: Laura Barreales y Sara Lévano
#Start date: 9th may 2025
#Version: 1.0

#Usage example
#symlinks.sh -d ruta/absoluta/del/directorio
version="version 1.0"
############################################################################################

# inicialización de logs vacíos
cat /dev/null > logs/stdout
cat /dev/null > logs/stderr

# explain the code
echo -e "Use: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" | tee -a logs/stdout

#Vamos a ver si se han proporcionado argumentos
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help"| tee -a logs/stdout
        exit 1
fi

#Manejar opciones y argumentos
while getopts ":d:hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script creates on symbolic link per file.fastq in the directory provided" | tee -a logs/stdout
                   echo -e "Usage: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" | tee -a logs/stdout
                   exit 0;;
                v) echo $version | tee -a logs/stdout
                   exit 0;;
                d) directory="$OPTARG"
                   echo "Procesando el directorio: $directory ..." | tee -a logs/stdout 
                    ;;
                :) echo "Option -$OPTARG requieres an argument" | tee -a logs/stdout #si pone solo -d no sirve
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" | tee -a logs/stdout
                   exit 1 ;;

        esac
done

for file in "$directory"/*; do
    echo -e "\n"
    sample=$(basename "$file") #Nos quedamos con el nombre de la muestra
    echo "This is $sample" 
    #Vamos a comprobar que el archivo no está vacío
    echo "Checking if the file exists and it's not empty" | tee -a logs/stdout
    if [ -s $file ]; then
            echo "File exist and is not empty"  | tee -a logs/stdout
    else
            echo "file is empty"  | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastq
    echo "Checking whether file extension is .fastq" | tee -a logs/stdout
    if [[ "$file" == *.fastq ]]; then
            echo "File extension is correct" | tee -a logs/stdout
    else
            echo "File extension not correct. Must be .fastq file" | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el archivo es legible y ejecutable
    echo "Checking $sample permissions" | tee -a logs/stdout

    if [[ -r "$file" && -x "$file" ]]; then
            echo "File is readable and executable" | tee -a logs/stdout
    elif [ -r "$file" ]; then
            echo "File is readable but not executable. Fixing." | tee -a logs/stdout
            chmod +x "$file"
            echo "fixed" | tee -a logs/stdout
    elif [ -x "$file" ]; then
            echo "File is executable but not readable" | tee -a logs/stdout
            chmod +r "$file"
            echo "Fieexd" | tee -a logs/stdout
    else
            echo "File is neither readable nor executable. Fixing" | tee -a logs/stdout
            chmod +rx "$file"
            echo "Fixed" | tee -a logs/stdout
    fi

    
    #Crearemos un link simbólico para cada archivo que se guardará en este directorio
    echo "Creating symbolic link for $sample ..." | tee -a logs/stdout
    ln -s $directory/$sample $sample   #Coge de la ruta del directorio dado cada archivo que haya en él y lo guarda como un symlink en este directorio
    mv $sample results/
    #Vamos a comprobar que haya ido bien la creación del link
    if [ $? != 0 ]; then
        echo "Something went wrong during symbolic links creation for file $sample" | tee -a logs/stderr
        exit 1
    fi
    echo "Symbolic link created for $file and named $sample" | tee -a logs/stdout
done
echo -e "\n"
echo "Symbolic links successfully created in $(pwd)" | tee -a logs/stdout

