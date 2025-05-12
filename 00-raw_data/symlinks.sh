# !/bin/bash
#En este script se crearan los links simbólicos para poder acceder a los archivos fasta de nuestas
#compañeras de clase.
#Author: Laura Barreales y Sara Lévano
#Start date: 9th may 2025
#Version: 1.0

#Usage example
#symlinks.sh -d ruta/absoluta/del/directorio
version="version 2.0"
############################################################################################

# inicialización de logs vacíos
cat /dev/null > logs/stdout
cat /dev/null > logs/stderr

# explain the code
{echo -e "Use: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" | tee -a logs/stdout

#Vamos a ver si se han proporcionado argumentos
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help"
        exit 1
fi

#Manejar opciones y argumentos
while getopts ":d:hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script creates on symbolic link per file.fastq in the directory provided" 
                   echo -e "Usage: \n\t-h for help\n\t-v for the version of the script \n\t-d ruta/absoluta/del/directorio to create de symbolic links" 
                   exit 0;;
                v) echo $version 
                   exit 0;;
                d) directory="$OPTARG"
                   echo "Procesando el directorio: $directory ..." 
                    ;;
                :) echo "Option -$OPTARG requieres an argument" >&2 #si pone solo -d no sirve
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" >&2
                   exit 1 ;;

        esac
done
} >> >(tee -a logs/all_samples.out) 2>> >(tee -a logs/all_samples.err)

{for file in "$directory"/*; do
    echo -e "\n"
    sample=$(basename "$file") #Nos quedamos con el nombre de la muestra
    name=${sample%.*} #quita la extensión del archivo, se usará para dar nombre a los logs de salida
    echo "This is $sample" 
    #Vamos a comprobar que el archivo no está vacío
    echo "Checking if the file exists and it's not empty" 
    if [ -s "$file" ]; then
            echo "File exist and is not empty"  
    else
            echo "file is empty" >&2
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastq
    echo "Checking whether file extension is .fastq" 
    if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" 
    else
            echo "File extension not correct. Must be .fastq file. Skipping $sample " >&2
            continue #Para que salte al siguiente archivo si no es un .fastq.gz
    fi 

    #Vamos a comprobar que el archivo es legible y ejecutable
    echo "Checking $sample permissions" 

    if [[ -r "$file" && -x "$file" ]]; then
            echo "File is readable and executable" 
    elif [ -r "$file" ]; then
            echo "File is readable but not executable. Fixing." >&2
            chmod +x "$file"
            echo "fixed" >&2
    elif [ -x "$file" ]; then
            echo "File is executable but not readable" >&2
            chmod +r "$file"
            echo "Fieexd" >&2
    else
            echo "File is neither readable nor executable. Fixing" >&2
            chmod +rx "$file"
            echo "Fixed" >&2
    fi

    
    #Crearemos un link simbólico para cada archivo que se guardará en este directorio
    echo "Creating symbolic link for $sample ..." 
    ln -s $file $sample   #Coge de la ruta del directorio dado cada archivo que haya en él y lo guarda como un symlink en este directorio
    mv $sample results/
    #Vamos a comprobar que haya ido bien la creación del link
    if [ $? != 0 ]; then
        echo "Something went wrong during symbolic links creation for file $sample" >&2
        exit 1  
    fi
    echo "Symbolic link created for $file and named $sample" 
done
} >> >(tee -a logs/${name}.out) 2>> >(tee -a logs/${name}.err)

echo -e "\n" >> >(tee -a logs/all_files.out)
echo "Symbolic links successfully created in $(pwd)" >> >(tee -a logs/all_files.out)

