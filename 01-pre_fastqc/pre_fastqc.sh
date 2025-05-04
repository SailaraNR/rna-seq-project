#Este script coge los archivos de la carpeta 00 qeu corresponden a las raw-data y analiza su calidad de secuencia
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Date: 4th May 2025
#Version: 1.0

#Usage example: pre-fastqc.sh
#No arguments expected, but -v and -h are available

#######################################################################

#Comprobar que el archivo es .fastqc

while getopts "hv" opt; do
        case $opt in
                h) echo "You asked for usage help"
                   echo "This script will use the check the quality of the previously downloaded raw data " 
                   echo -e "You won't need to add any flags, but -v and -h are available to ask for script version \n or help" ;;
                v) echo "Version 1.0" 
                   exit 1 ;;
                ?) echo "Invalid option or missing argument" >&2
                   echo "Usage example: $0 -u URL" >&2
                   exit 1;;
        esac
done

for files in ../00-raw_data; do
    #Vamos a comprobar que el archivo es "readable" y ejecutable
    echo "Cheking $file permissions"
    if [ -r $file ]; then
        if [ -x $file ]; then
                echo "This file is readable and executable"
        else
                echo "This file is readable but not executable"
                chmod +x $file
                echo "Now is ready"
    elif [ -x $file ]; then
        echo "This file is not readable but executable"
        chmod +r $file
        echo "Now is ready"
    else
        echo "This file es not readable not executable"
        chmod +rx $file
        echo "Now is ready"

    fi
    #También se comprobará si tiene la extensión correcta (.fastqc)
    echo "Checking whether file extension is .fastq"
    if [ "$file" == *.fastq ]; then
        echo "File extension is correct"
    else
        echo "File extension not correct. Must be .fastq file" >&2
        exit 1
    fi
    #También se comprobará si el archivo está vacío
    echo "Cheking if the file exists and it's not empty"
    if [ -s $file ]; then
        echo "File exist and is not empty" 
    else
        echo "$file is empty"  >&2
        exit 1
    fi 
    #Ahora vamos a comprobar que la calidad de la secuencia con fastqc
    echo "Cheking quality score in $file"
    fastqc $file -o ./results/fastqc_$file.fastqc #-o guarda los archivos en esa ruta, que en este caso es la carpeta results de la carpeta 01
    echo "Fastqc analysis completed for $file"
done 

echo "Now Multiqc will run in order to make a summary of the analysis"
multiqc ./results -n multiqc_analysis.html #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc
# se habrá creado un archivo html con un resumen de el análisis de fastqc
# -n es para ponerle nombre al archivo
echo -e "Multiqc has finished its job. \n You can now check for an .html file in $pwd \n to find the analysis"

######################################################
#Voy a dejar este enlace aquí, ya lo quitaremos.
#Es para usar multiqc https://docs.seqera.io/multiqc/getting_started/running_multiqc