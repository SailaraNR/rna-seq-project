#Este script acepta una URL que corresponde a un archivo con una lista de las accesion de las muestras que se queiren descargar
#Se puede obetener este archivo desde la web de SRA correspondiente al experimento
#Author: Laura Barreales y Sara Lévano
#Date: 1st may 2025
#Version: 1.1

#Usage example
#00-data.sh -u URL

#ARguments expected: -u URL

############################################################################################

#comprobar que se ha asignado bien la URL y asignarla a una variable
echo "Let's check if the argument is correctly given"
while getopts "u:hv" opt; do
        case $opt in
                h) echo "You asked for usage help" 
                   echo "This script downloads rawdata if a file.txt with SRR accessions is given"
                   echo "You will need to execute 00-data.sh -u file.txt" 
                   echo "Then the raw sequences you want will be downloaded";;
                v) echo "Version 1.1" ;;
                u) URL="$OPTARG" ;;
                ?) echo "Invalid option or missing argument" >&2
                   echo "Usage example: $0 -u URL" >&2\
                   exit 1;;
        esac
done


#Vamos a descargar el archivo de la URL dada
echo "Downloading the file"
curl -L -o file "$URL"

#VAmos a comprobar que se ha descargado bien
echo "Checking if the file has been succesfully downloaded"
if [ $? -ne 0 ]; then
        echo "Downloading error" >&2
        exit 1
fi

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


#Vamos a comprobar que el archivo no está vacío
echo "Cheking if the file exists and it's not empty"
if [ -s $file ]; then
        echo "File exist and is not empty" 
else
        echo "Downloaded file is empty"  >&2
        exit 1
fi 


#Vamos a comprobar que el arhivo es .txt
echo "Checling whether file extension is .txt"
if [ "$file" == *.txt ]; then
        echo "File extension is correct"
else
        echo "File extension not correct. Must be .txt file" >&2
        exit 1
fi

#Se leerán las accesion del archivo una a una y se descargaran las raw data correspondientes
echo "Processing accesions from $FILENAME"
while read -r ACCESSION; do
  [[ -z "$ACCESSION" ]] && continue # comprueba que la linea que lee está o no vacía
  echo "Downloading $ACCESSION"
  prefetch "$ACCESSION" #descarga .sra files
  fasterq-dump "$ACCESSION" --split-files #convierte .sra files en .fastq y los separa en dos archivos si se trata de lecturas pareadas
done < "$FILENAME"

echo "Finished. FASTQ files downloaded in $(pwd)"
ls #para ver los archivos
