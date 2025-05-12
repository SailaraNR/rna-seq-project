# !/bin/bash

#en este script voy a modificar las redirecciones, pero como el otro funciona creo uno aparte

#Este script coge los archivos de la carpeta 00 qeu corresponden a las raw-data y analiza su calidad de secuencia
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 4th May 2025
version="Version 2.0 "

#Usage example: pre-fastqc.sh
#No arguments expected, but -v and -h are available
#link al manual de FastQC: https://github.com/s-andrews/FastQC/blob/master/fastqc
#link al manual de MultiQC. No es de github, pero he accedido desde github:https://docs.seqera.io/multiqc/getting_started/running_multiqc 

##These must be the input paht and output path in order to correctely running the script
input_path="../00-raw_data/results"
output_path="./results"
#######################################################################

# inicialización de logs vacíos
cat /dev/null > logs2/*.out
cat /dev/null > logs2/*.err

#explain the code
echo -e "Use -h for help or -v for the version. \nUse -i to introduce input directory and -o to introduce output directory for fastqc outputfiles" | tee -a logs/help.out
# getops y while para la ayuda, la versión y el ejemplo de uso

#Vamos a ver si se han proporcionado el número de argumentos correctos (0 o 1)
if [ $# -gt 1 ]; then
        echo -e "More than 1 option provided. \nRun $0 -h for usage help"| tee -a logs/help.err
        exit 1
fi

#Manjar opciones y argumentos
while getopts ":hvo:i:" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script checks the quality of the previously downloaded raw sequences" | tee -a logs/help.out
            echo -e "Usage example:\n\t./$0 -h: Provides help \n\t./$0 -v: Tells script's version\n\t./$0: checks quality score of raw sequences" \
            "\n\t./$0 -i /path/to/input/direcyoty \n\t./$0 -o /path/to/output/direcyoty" | tee -a logs/help.out
            exit 0;;
        v) echo "$version" | tee -a logs/help.out
            exit 0;;
        i) input_path="$OPTARG"
            echo "Procesando el directorio: $input_path ..." ;;
        o) output_path="$OPTARG"
            echo "Procesando el directorio: $output_path ..." ;;
        :) echo "Option -$OPTARG requieres an argument" >&2 #si pone solo -d no sirve
            exit 1 ;;    
        \?) echo -e "Invalid option\n\ 
            Usage example: $0 " | tee -a logs/help.err
            exit 1;; #si pone simbolos raros no sirve
    esac
done


# Comprobar que los directorios existan
echo -e "\nChecking if output directory exists" | tee -a logs/help.out
if [ ! -d $input_path ]
then 
     echo "Input path exists"
else 
     if [ "$(ls $input_path)" ]
     then  
         echo "Directory is not empty"
     else 
         echo "Directory is empty"
     fi
fi

echo -e "\nChecking if output directory exists" | tee -a logs/help.out
if [ ! -d $output_path ]
then 
     echo "Input path exists"
else 
     if [ "$(ls $output_path)" ]
     then  
         echo "Directory is not empty"
     else 
         echo "Directory is empty"
     fi
fi



{ for file in "$input_path"/*; do
        #Vamos a comprobar que el archivo no está vacío
        sample=$(basename "$file") #Nos quedamos con el nombre de la muestra para los echo
        name=${sample%.*} #quita la extensión del archivo, se usará para dar nombre a los files de salida
        echo -e "\nChecking if $sample exists and it's not empty" 
        if [[ -f "$file" && -s "$file" ]]; then
            echo "$sample exist and is not empty"  
         else
            echo "$sample is not a file or is empty"  >&2
            exit 1
        fi 

        #Vamos a comprobar que el arhivo es .fastqc
        echo "Checking whether file extension is .fastq.gz" 
        if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" 
            #Vamos a comprobar que el archivo es legible. Solo hace falta que sea legible
            echo "Checking $sample permissions" 
            if [[ ! -r "$file" ]]; then
                echo "File is not readable. Fixing..." >&2
                chmod +r "$file" && echo "Permissions fixed" || echo "Could not change permisions" >&2
                    continue
            fi
            
            #Ahora vamos a comprobar que la calidad de la secuencia con fastqc
            echo "Checking quality score in $sample" 
            #creo que podría ponerse así 
            fastqc $file -o $output_path   #porque $file debería ser la ruta absoluta del archivo (al menos así era en el archivo symlinks)
           
            #quitaría la linea de abajo, está mal para este bucle concreto
            #fastqc $input_path/*.fastq -o $output_path | tee -a logs/${sample}.out #-o guarda los archivos en esa ruta, que en este caso es la subcarpeta results
            if  [ $? -eq 0 ]; then
                echo "Fastqc analysis completed for $sample" 
            else
                echo "Fastqc analysis failed for $sample" >&2
            fi    
        else
            echo "this is not .fastq.gz, moving to the next..." >&2
            continue
            #Si funciona lo de guardar los symlinks en 00/results no hace falta esto, pero se podría dejar como cortafuegos
        fi;
done } >> >(tee -a logs/${name}.out) 2>> >(tee -a logs/${name}.err)

{ echo -e "\nNow Multiqc will run in order to make a summary of the analysis" 
echo "But first it will check if files' extensions are .html after fastqc analysis" 
#Vamos a comprobar que el arhivo es .html
#for file in "$output_path"/*; do
#    echo "Checking whether file extension is .html" 
#    if [[ "$file" == *.html ]]; then
#            echo "File extension is correct" 
#    else
#            echo "File extension not correct. Must be .html file" 
#            exit 1
#    fi 
#done 
#multiqc ./results -n multiqc_analysis2.html | tee -a logs/${sample}.out \ 2> >(tee -a logs/${sample}.err >&2) #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc

#A lo mejor así queda más elegante
multiqc $output_path -n multiqc_analysis2.html #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc
# se habrá creado un archivo html con un resumen de el análisis de fastqc
# -n es para ponerle nombre al archivo
echo -e "Multiqc has finished its job. \nYou can now check for an .html file in $(pwd) \nto find the analysis" 
} >> >(tee -a logs/multiqc.out) 2>> >(tee -a logs/multiqc.err)

######################################################
#Voy a dejar este enlace aquí, ya lo quitaremos.
#Es para usar multiqc https://docs.seqera.io/multiqc/getting_started/running_multiqc
