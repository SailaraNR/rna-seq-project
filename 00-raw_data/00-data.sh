# !/bin/bash
#Este script acepta una URL que corresponde a un archivo con una lista de las accesion de las muestras que se queiren descargar
#Se puede obetener este archivo desde la web de SRA correspondiente al experimento
#Author: Laura Barreales y Sara Lévano
#Start date: 1st may 2025
#Version: 2.0.1

#Usage example
#00-data.sh -f SRR_Acc_List.txt

#Arguments expected: -f file.txt

############################################################################################

# inicialización de logs vacíos
cat /dev/null > logs/stdout
cat /dev/null > logs/stderr

# explain the code
echo "Use -f for the file, -h for help and -v for the version" | tee -a logs/stdout
# getops y while para la ayuda, la versión y el ejemplo de uso

#Vamos a ver si se han proporcionado argumetnos
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help"| tee -a logs/stdout
        exit 1
fi

#Manejar opciones y argumentos
while getopts ":f:hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script downloads rawdata if a file.txt with SRR accessions is given" | tee -a logs/stdout
                   echo -e "Usage example:\n\t sh $0 -h: Provides help \n\t sh $0 -v: Tells script's version\n\t sh $0 -f file.txt: Downloads rawdata from de accessions in file.txt" | tee -a logs/stdout
                   exit 0;;
                v) echo "Version 1.0" | tee -a logs/stdout
                   exit 0;;
                f) echo "Abriendo el archivo -f $file ..." | tee -a logs/stdout 
                   file="$OPTARG" ;;
                :) echo "Option -$OPTARG requieres an argument" | tee -a logs/stdout #si pone solo -f no sirve
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" | tee -a logs/stdout
                   exit 1 ;;

        esac
done

#Vamos a comprobar que el archivo no está vacío
echo "Checking if the file exists and it's not empty" | tee -a logs/stdout
if [ -s $file ]; then
        echo "File exist and is not empty"  | tee -a logs/stdout
else
        echo "file is empty"  | tee -a logs/stderr
        exit 1
fi 

#Vamos a comprobar que el arhivo es .txt
echo "Checking whether file extension is .txt" | tee -a logs/stdout
if [ "$file" == *.txt ]; then
        echo "File extension is correct" | tee -a logs/stdout
else
        echo "File extension not correct. Must be .txt file" | tee -a logs/stderr
        exit 1
fi 

#Vamos a comprobar que el archivo es legible y ejecutable
echo "Checking $file permissions" | tee -a logs/stdout

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

#Se leerán las accesion del archivo una a una y se descargaran las raw data correspondientes
echo "Processing accesions from $file..." | tee -a logs/stdout
while read -r ACCESSION; do
  [[ -z "$ACCESSION" ]] && continue # comprueba que la linea que lee está o no vacía
  echo "Downloading $ACCESSION" | tee -a logs/stdout
  prefetch "$ACCESSION" #descarga .sra files
  fasterq-dump --split-files "$ACCESSION" #convierte .sra files en .fastq y los separa en dos archivos si se trata de lecturas pareadas
done < "$file"

echo "Finished. FASTQ files downloaded in $(pwd)" | tee -a logs/stdout
ls #para ver los archivos
