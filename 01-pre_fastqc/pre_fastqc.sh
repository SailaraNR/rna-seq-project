#Este script coge los archivos de la carpeta 00 qeu corresponden a las raw-data y analiza su calidad de secuencia
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 4th May 2025
#Version: 1.2

#Usage example: pre-fastqc.sh
#No arguments expected, but -v and -h are available

#######################################################################

#explain the code
echo -e "Use -h for help or -v for the version. \nBut do not introduce any arguments or options for checking the quality score" | tee -a logs/stdout
# getops y while para la ayuda, la versión y el ejemplo de uso

#Vamos a ver si se han proporcionado el número de argumentos correctos (0 o 1)
if [ $# -gt 1 ]; then
        echo -e "More than 1 option provided. \nRun $0 -h for usage help"| tee -a logs/stdout
        exit 1
fi

#Manjar opciones y argumentos
while getopts "hv" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script checks the quality of the previously downloaded raw sequences" | tee -a logs/stdout
            echo -e "Usage example:\n\t sh $0 -h: Provides help \n\t sh $0 -v: Tells script's version\n\tsh $0: checks quality score of raw sequences" | tee -a logs/stdout
            exit 0;;
        v) echo "Version 1.2" | tee -a logs/stdout
            exit 0;;
        \?) echo -e "Invalid option\n\ 
            Usage example: $0 " | tee -a logs/stdout
            exit 1;; #si pone simbolos raros no sirve
    esac
done

for file in ../00-raw_data/results; do
    #Vamos a comprobar que el archivo no está vacío
    echo "Cheking if the file exists and it's not empty" | tee -a logs/stdout
    if [ -s $file ]; then
            echo "File exist and is not empty"  | tee -a logs/stdout
    else
            echo "file is empty"  | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastqc
    echo "Checling whether file extension is .fastq" | tee -a logs/stdout
    if [ "$file" == *.fastq ]; then
            echo "File extension is correct" | tee -a logs/stdout
    else
            echo "File extension not correct. Must be .fastq file" | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el archivo es legible y ejecutable
    echo "Cheking $file permissions" | tee -a logs/stdout

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
    #Ahora vamos a comprobar que la calidad de la secuencia con fastqc
    echo "Cheking quality score in $file"
    fastqc $file -o ./results/fastqc_$file.html #-o guarda los archivos en esa ruta, que en este caso es la carpeta results de la carpeta 01
    echo "Fastqc analysis completed for $file"
done 

echo "Now Multiqc will run in order to make a summary of the analysis"
echo "Bur first it will check if files' extensions are .html after fastqc analysis"
#Vamos a comprobar que el arhivo es .html
    echo "Checling whether file extension is .html" | tee -a logs/stdout
    if [ "$file" == *.html ]; then
            echo "File extension is correct" | tee -a logs/stdout
    else
            echo "File extension not correct. Must be .fastq file" | tee -a logs/stderr
            exit 1
    fi 
multiqc ./results -n multiqc_analysis.html #multiqc se va a ejecutar en el 01/results que es donde se han guardado las secuencias fastqc
# se habrá creado un archivo html con un resumen de el análisis de fastqc
# -n es para ponerle nombre al archivo
echo -e "Multiqc has finished its job. \n You can now check for an .html file in $pwd \n to find the analysis"

######################################################
#Voy a dejar este enlace aquí, ya lo quitaremos.
#Es para usar multiqc https://docs.seqera.io/multiqc/getting_started/running_multiqc