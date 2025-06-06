# !/bin/bash
#Este script coge los archivos de la carpeta 00 qeu corresponden a las raw-data y analiza su calidad de secuencia
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 4th May 2025
version="Version 2.0"

#Usage example: pre-fastqc.sh
#No arguments expected, but -v and -h are available
#link al manual de FastQC: https://github.com/s-andrews/FastQC/blob/master/fastqc
#link al manual de MultiQC. No es de github, pero he accedido desde github:https://docs.seqera.io/multiqc/getting_started/running_multiqc 

#######################################################################

# inicialización de logs vacíos
cat /dev/null > logs/stdout
cat /dev/null > logs/stderr

#explain the code
echo -e "Use -h for help or -v for the version. \nBut do not introduce any arguments or options for checking the quality score" | tee -a logs/stdout
# getops y while para la ayuda, la versión y el ejemplo de uso

#Vamos a ver si se han proporcionado el número de argumentos correctos (0 o 1)
if [ $# -gt 1 ]; then
        echo -e "More than 1 option provided. \nRun $0 -h for usage help"| tee -a logs/stderr
        exit 1
fi

#Manjar opciones y argumentos
while getopts "hv" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script checks the quality of the previously downloaded raw sequences" | tee -a logs/stdout
            echo -e "Usage example:\n\tsh $0 -h: Provides help \n\tsh $0 -v: Tells script's version\n\tsh $0: checks quality score of raw sequences" | tee -a logs/stdout
            exit 0;;
        v) echo "$version" | tee -a logs/stdout
            exit 0;;
        \?) echo -e "Invalid option\n\ 
            Usage example: $0 " | tee -a logs/stderr
            exit 1;; #si pone simbolos raros no sirve
    esac
done

input_path="../00-raw_data/results"
output_path="/results/"

# Check if the target is not a directory
echo -e "\nChecking if input directory exists"
if [ ! -d "$input_path" ]; then
        echo "This is not a directory" | 2> logs/stderr
        exit 1
fi

{ for file in "$input_path"/*; do
        #Vamos a comprobar que el archivo no está vacío
        sample=$(basename "$file") #Nos quedamos con el nombre de la muestra
        echo -e "\nChecking if $sample exists and it's not empty" | tee -a logs/stdout
        if [[ -f "$file" && -s "$file" ]]; then
            echo "$sample exist and is not empty"  | tee -a logs/stdout
         else
            echo "$sample is not a file or is empty"  | tee -a logs/stderr
            exit 1
        fi 

        #Vamos a comprobar que el arhivo es .fastqc
        echo "Checking whether file extension is .fastq.gz" | tee -a logs/stdout
        if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" | tee -a logs/stdout
            #Vamos a comprobar que el archivo es legible. Solo hace falta que sea legible
            echo "Checking $sample permissions" | tee -a logs/stdout
            if [[ ! -r "$file" ]]; then
                echo "File is not readable. Fixing..." | tee -a logs/stdout
                chmod +r "$file" && echo "Permissions fixed" | tee -a logs/stdout || echo "Could not change permisions" | tee -a logs/stderr 
                    continue
            fi
            
            #Ahora vamos a comprobar que la calidad de la secuencia con fastqc
            echo "Checking quality score in $sample" | tee -a logs/stdout
            #creo que podría ponerse así 
            fastqc $file -o $output_path 2>> logs/stderr | tee -a logs/stdout #porque $file debería ser la ruta absoluta del archivo (al menos así era en el archivo symlinks)
           
            #quitaría la linea de abajo, está mal para este bucle concreto
            #fastqc $input_path/*.fastq -o $output_path | tee -a logs/stdout #-o guarda los archivos en esa ruta, que en este caso es la subcarpeta results
            if  [ $? -eq 0 ]; then
                echo "Fastqc analysis completed for $sample" | tee -a logs/stdout
            else
                echo "Fastqc analysis failed for $sample"
            fi    
        else
            echo "this is not .fastq.gz, moving to the next..." 
            continue
            #Si funciona lo de guardar los symlinks en 00/results no hace falta esto, pero se podría dejar como cortafuegos
        fi;
done } 2>> logs/stderr | tee -a logs/stdout

echo -e "\nNow Multiqc will run in order to make a summary of the analysis" | tee -a logs/stdout
echo "But first it will check if files' extensions are .html after fastqc analysis" | tee -a logs/stdout
#Vamos a comprobar que el arhivo es .html
#for file in "$output_path"/*; do
#    echo "Checking whether file extension is .html" | tee -a logs/stdout
#    if [[ "$file" == *.html ]]; then
#            echo "File extension is correct" | tee -a logs/stdout
#    else
#            echo "File extension not correct. Must be .html file" | tee -a logs/stderr
#            exit 1
#    fi 
#done 
#multiqc ./results -n multiqc_analysis.html #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc

#A lo mejor así queda más elegante
multiqc $output_path -n multiqc_analysis.html 2>> logs/stderr | tee -a logs/stdout #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc
# se habrá creado un archivo html con un resumen de el análisis de fastqc
# -n es para ponerle nombre al archivo
echo -e "Multiqc has finished its job. \nYou can now check for an .html file in $(pwd) \nto find the analysis" | tee -a logs/stdout

######################################################
#Voy a dejar este enlace aquí, ya lo quitaremos.
#Es para usar multiqc https://docs.seqera.io/multiqc/getting_started/running_multiqc
