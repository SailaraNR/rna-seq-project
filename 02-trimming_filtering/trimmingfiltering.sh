# !/bin/bash
#Este script coge los archivos con las raw data, las analiza, trima y elimina secuencias
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 6th May 2025
version="Version: 2.1"

#Usage example: trimmingfiltering.sh
#No arguments expected, but -v and -h are available

##The directories must be:

#input_dir="../00-raw_data/results"
#output_dir="./results"
#fastqc_dir="./results/fastqc"
#multiqc_dir="./results/multiqc"

#TENEMOS QUE CREAR UNA CARPETA PARA MULTIQC Y FASTQC DENTRO DE LOS LOGS

#######################################################################

#explain the code
{ echo -e "Use -h for help or -v for the version. \nBut do not introduce any arguments or options for for trimming and filtering sequences" | tee -a logs/stdout

#Parámetros por defecto de Trimmomatic por si el usuario no pone ninguno
#Si el usuario pone alguno, estas variables se reescribiran gracias al getopts 
SLIDINGWINDOW="4:20"
MINLEN="36"
LEADING="10"
TRAILING="10"

#Manejar opciones y argumentos
while getopts "hvm:w:l:t:i:o:f:c:" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script trims and filters raw sequences" 
            echo -e "Usage example:\n\t./$0 -h: Provides help \n\t./$0 -v: Tells script's version\n\t./$0 -w -m -l -t -i -o -f -c: trims and filters raw sequences" 
            exit 0;;
        v) echo "$version" 
            exit 0;;
        w) SLIDINGWINDOW="$OPTARG" ;;
        m) MINLEN="$OPTARG" ;;
        l) LEADING="$OPTARG" ;;
        t) TRAILING="$OPTARG" ;;
        i) input_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        f) fastqc_dir="$OPTARG" ;;
        c) multiqc_dir="$OPTARG" ;;
        \?) echo -e "Invalid option\n"\ 
            "Usage example: $0 " >&2
            exit 1;; #si pone simbolos raros no sirve
    esac
done
echo "Usando SLIDINGWINDOW=$SLIDINGWINDOW y MINLEN=$MINLEN" 
} >> >(tee -a logs/all_files.out) 2>> >(tee -a logs/all_files.err)

for file in "$input_dir"/*; do #Los archivos fastq para trimar y filtrar se encuentran en esa carpeta. $file debería ser una ruta ansoluta al file
    sample=$(basename "$file") #Nos quedamos con el nombre de la muestra. Falta quitarle el .fastq, ya veremos como hacerlo en clase
    name=${sample%.*} #quitar la extension del archivo para normbrar los logs
    #Vamos a comprobar que el archivo no está vacío
    { echo "Cheking if the "$sample" exists and it's not empty" 
    if [ -s "$file" ]; then
            echo "File exist and is not empty"  
    else
            echo "file is empty"  >&2
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastq
    echo "Checking whether file extension is .fastq" 
    if [[ "$file" == *.fastq ]]; then
            echo "File extension is correct" 
    else
            echo "File extension not correct. Must be .txt file" >&2
            exit 1
    fi 

    #Vamos a comprobar que el archivo es legible y ejecutable
    echo "Checking "$sample" permissions" 

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
    } >> >(tee -a logs/${name}.out) 2>> >(tee -a logs/${name}.err)
    #Ahora vamos a trimar. 
    #Necesitamos saber como se llaman los archivos fastq y cómo se diferencian entre pares de secuencias
    #Va a sacar 4 archivos, dos de ellos pasarán la selección tras el trimado (paired) y otros dos no (unpaired). Esto ocurre porque son secuencias pareadas
    echo "Processing "$sample"..."
    name=${sample%.*} #Para quitar la extensión y poder poner el nombre de la muestra en los archivos input y putput
    {
    trimmomatic PE \
      ${name}_R1.fastq ${name}_R2.fastq \
      ${name}_R1_paired.fastq ${name}_R1_unpaired.fastq \
      ${name}_R2_paired.fastq ${name}_R2_unpaired.fastq \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:$SLIDINGWINDOW \
      LEADING:$LEADING TRAILING:$TRAILING \
      MINLEN:$MINLEN
     } >> >(tee -a logs/${name}_TRIM.out) 2>> >(tee -a logs/${name}_TRIM.err) 
      #MINLEN: Es la cantidad mínima de bases que pedimos tras el trimado, yo pondría más
      #LEADING: Es la calidad mínima que debe tener la secuenciación de cada base del principio, si no llega se corta
      #TRAILING: Es la calidad mínima que debe tener la secuenciación de cada base del final, si no llega se corta
      #SLIDINGWINDOW:<windowSize>:<requiredQuality>  Mira en la ventana de bases cual es la calidad media, si no supera un mínimo corta. Compatible con TRAILING y LEADING?
      #ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>  Preguntar si esas secuencias tienen adaptadores
      #TOPHRED33/64 This (re)encodes the quality part of the FASTQ file to base 33 or 64.
      #NO utilizaría MAXINFO 
      
      #Hacaer Fastqc para cada muestra
      {
      echo "Executing FastQC for paired files of $sample" | tee -a logs/stdout
      fastqc "$output_dir/${name}_R1_paired.fastq" "$output/${name}_R2_paired.fastq" -o "$fastqc_dir"
      } >> >(tee -a logs/fastqc/${name}_fastqc.out) 2>> >(tee -a logs/fastqc/${name}_fastqc.err)
done


# Hacer multiqc para integrar el análisis del fastqc
{ echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$multiqc_dir" -n "multiqc_trimmed_analysis.html" -o "$multiqc_dir" } >> >(tee -a logs/multiqc/multiqc.out) 2>> >(tee -a logs/multiqc/multiqc.err)
   

echo -e "Analysis completed." >> >(tee -a logs/all_files.out)
