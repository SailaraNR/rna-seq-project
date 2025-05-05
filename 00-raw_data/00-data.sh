#Este script acepta una URL que corresponde a un archivo con una lista de las accesion de las muestras que se queiren descargar
#Se puede obetener este archivo desde la web de SRA correspondiente al experimento
#Author: Laura Barreales y Sara Lévano
#Date: 1st may 2025
#Version: 1.3

#Usage example
#00-data.sh -f SRR_Acc_List.txt

#Arguments expected: -f file.txt

############################################################################################

#explain the code
echo "Use -f for the file, -h for help and -v for the version"
# getops y while para la ayuda, la versión y el ejemplo de uso
{
while getopts "hv" opt; do
        case $opt in
                h) echo -e "You asked for usage help\n"\
                   "This script downloads rawdata if a file.txt with SRR accessions is given" | tee -a logs/stdout
                   exit 0;;
                v) echo "Version 1.3" | tee -a logs/stdout
                   exit 0;;
                f) echo "Abriendo el archivo $2 ..." | tee -a logs/stdout
                   exit 0;;
                ?) echo -e "Invalid option or missing argument\n"\
                   "Usage example: $0 SRR_Acc_List" | tee -a logs/stderr
                   exit 1;;
        esac
done
} >> logs/stdout 2>> logs/stderr

# Mostrar uso correcto del script
if [ -z "$2" ]; then
    echo "Error: No se ha proporcionado un archivo.\n"\
    "Uso: $0 -f SRR_Acc_List.txt" | tee -a logs/stderr
    exit 1
fi >> logs/stdout 2>> logs/stderr

#Vamos a comprobar que el archivo es legible y ejecutable
echo "Cheking $file permissions"

{
if [ -r $file ]; then
        if [ -x $file ]; then
                echo "This file is readable and executable" | tee -a logs/stdout
        else
                echo "This file is readable but not executable" | tee -a logs/stdout
                chmod +x $file
                echo "Now is ready" | tee -a logs/stdout
elif [ -x $file ]; then
        echo "This file is not readable but executable" | tee -a logs/stdout
        chmod +r $file
        echo "Now is ready" | tee -a logs/stdout
else
        echo "This file es not readable not executable" | tee -a logs/stdout
        chmod +rx $file
        echo "Now is ready" | tee -a logs/stdout

fi
} >> logs/stdout 2>> logs/stderr

#Vamos a comprobar que el archivo no está vacío
{
echo "Cheking if the file exists and it's not empty"
if [ -s $file ]; then
        echo "File exist and is not empty"  | tee -a logs/stdout
else
        echo "file is empty"  | tee -a logs/stderr
        exit 1
fi 
} >> logs/stdout 2>> logs/stderr

#Vamos a comprobar que el arhivo es .txt
{
echo "Checling whether file extension is .txt"
if [ "$file" == *.txt ]; then
        echo "File extension is correct" | tee -a logs/stdout
else
        echo "File extension not correct. Must be .txt file" | tee -a logs/stderr
        exit 1
fi 
} >> logs/stdout 2>> logs/stderr

#Se leerán las accesion del archivo una a una y se descargaran las raw data correspondientes
{
echo "Processing accesions from $file..."
while read -r ACCESSION; do
  [[ -z "$ACCESSION" ]] && continue # comprueba que la linea que lee está o no vacía
  echo "Downloading $ACCESSION" | tee -a logs/stdout
  prefetch "$ACCESSION" #descarga .sra files
  fasterq-dump "$ACCESSION" --split-files #convierte .sra files en .fastq y los separa en dos archivos si se trata de lecturas pareadas
done < "$file"
} >> logs/stdout 2>> logs/stderr

echo "Finished. FASTQ files downloaded in $(pwd)" | tee -a logs/stdout
ls #para ver los archivos
