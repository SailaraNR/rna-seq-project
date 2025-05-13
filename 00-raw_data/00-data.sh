# !/bin/bash
#Este script acepta una URL que corresponde a un archivo con una lista de las accesion de las muestras que se queiren descargar
#Se puede obetener este archivo desde la web de SRA correspondiente al experimento
#Author: Laura Barreales y Sara Lévano
#Start date: 1st may 2025
version= Version 2.3

#Usage example
#00-data.sh -f SRR_Acc_List.txt

#Arguments expected: -f file.txt

############################################################################################

# inicialización de logs vacíos
cat /dev/null > logs/*.out
cat /dev/null > logs/*.err

# explain the code
{ echo "Use -f for the file, -h for help and -v for the version" 

#Vamos a ver si se han proporcionado argumetnos
if [ $# -eq 0 ]; then
        echo -e "No options provided. \nUse $0 -h for help" >&2
        exit 1
fi

#Manejar opciones y argumentos
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
                :) echo "Option -$OPTARG requieres an argument" >&2 #si pone solo -f no sirve
                   exit 1 ;;
                \?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List.txt" >&2
                   exit 1 ;;

        esac
done
} 2>> >(tee -a logs/all_samples.err) >> >(tee -a logs/all_samples.out)

#Vamos a comprobar que el archivo no está vacío
{ echo "Checking if the file exists and it's not empty" 
if [ -s $file ]; then
        echo "File exist and is not empty"  
else
        echo "file is empty"  >&2
        exit 1
fi 

#Vamos a comprobar que el arhivo es .txt
echo "Checking whether file extension is .txt" 
if [ "$file" == *.txt ]; then
        echo "File extension is correct" 
else
        echo "File extension not correct. Must be .txt file" >&2
        exit 1
fi 

#Vamos a comprobar que el archivo es legible y ejecutable
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
        echo "Fieexd" >&2
else
        echo "File is neither readable nor executable. Fixing" >&2
        chmod +rx "$file"
        echo "Fixed" >&2
fi 
} 2>> >(tee -a logs/all_samples.err) >> >(tee -a logs/all_samples.out)

#Se leerán las accesion del archivo una a una y se descargaran las raw data correspondientes
{ echo "Processing accesions from $file..." 
while read -r ACCESSION; do
  name=$(sample%.*)
  [[ -z "$ACCESSION" ]] && continue # comprueba que la linea que lee está o no vacía
  echo "Downloading $ACCESSION" 
  prefetch "$ACCESSION" #descarga .sra files
  fasterq-dump --split-files "$ACCESSION" #convierte .sra files en .fastq y los separa en dos archivos si se trata de lecturas pareadas
done < "$file"
rm -r *.sra #para borrar los sra files descargados con prefetch que no sirven para nada

echo "Finished. FASTQ files downloaded in $(pwd)" 
} 2>> >(tee -a logs/all_${name}.err) >> >(tee -a logs/all_${name}.out)
