# !/bin/bash
#Este script coge los archivos con las raw data, las analiza, trima y elimina secuencias
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 6th May 2025
version="Version: 3.0"

#Usage example: trimmingfiltering.sh
#No arguments expected, but -v and -h are available

##The directories must be:

#input_dir="../00-raw_data/results"
#output_dir="./results"
#fastqc_dir="./results/fastqc"
#multiqc_dir="./results/multiqc"


#TENEMOS QUE CREAR UNA CARPETA PARA MULTIQC Y FASTQC DENTRO DE LOS LOGS

#######################################################################
# inicialización de logs vacíos
cat /dev/null > logs/*

#explain the codec
{ echo -e "Use -h for help or -v for the version."

#Parámetros por defecto de Trimmomatic por si el usuario no pone ninguno
#Si el usuario pone alguno, estas variables se reescribiran gracias al getopts 
SLIDINGWINDOW="4:20"
MINLEN="40" #ente 36-50 suele estar bien
LEADING="3" #como elimina desde el principio donde la secuenciación tiene baja calidad, evita cortar demaisado
TRAILING="20" #elimina desde atrás pb con phredscore menor a 20
ILLUMINACLIP="TruSeq3-PE-2.fa:2:30:10" #contiene el archivo de los adaptadores: el máximo de mismatch que se aceptan:palindromeClipthreshold:simpleCLipThreshold

#Manejar opciones y argumentos
while getopts "hvm:w:l:t:i:o:f:c:p:" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script trims and filters raw sequences" 
            echo -e "Usage example:\n\t$0 -h: Provides help \n\t$0 -v: Tells script's version\n\t$0 -w -m -l -t -i -o -f -c: trims and filters raw sequences" 
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
        p) ILLUMINACLIP="$OPTARG" ;;
        \?) echo -e "Invalid option\nUsage example: $0 " >&2
            exit 1;; #si pone simbolos raros no sirve
    esac
done
echo "Usando SLIDINGWINDOW=$SLIDINGWINDOW y MINLEN=$MINLEN" 
} 2>> >(tee -a logs/all_files.err) >> >(tee -a logs/all_files.out) 
#mkdir -p "$output_dir" "$fastqc_dir" "$multiqc_dir" 


#Comprobar que los directorios existan
echo -e "\nChecking if output directory exists" 
if [ ! -d "$input_dir" ]
then 
     echo "$input_dir does not path exist"
fi
#Comprobar que los directorios existan
echo -e "\nChecking if output directory exists" 
if [ ! -d "$output_dir" ]
then 
     echo "$output_dir does not path exist"
fi
#Comprobar que los directorios existan
echo -e "\nChecking if output directory exists" 
if [ ! -d "$fastqc_dir" ]
then 
     echo "$fastqc_dir does not path exist"
fi
#Comprobar que los directorios existan
echo -e "\nChecking if output directory exists" 
if [ ! -d "$multiqc_dir" ]
then 
     echo "$multiqc_dir does not path exist"
fi

#crea los directorios que especifíca
for file in "$input_dir"/*_1.fastq.gz; do #Los archivos fastq para trimar y filtrar se encuentran en esa carpeta. $file debería ser una ruta ansoluta al file
    sample=$(basename "$file" "_1.fastq.gz") #Nos quedamos con el nombre de la muestra.  
    #name=${sample%.*} #Para quitar la extensión y poder poner el nombre de la muestra en los archivos input y putput
    file1="${input_dir}/${sample}_1.fastq.gz"
    file2="${input_dir}/${sample}_2.fastq.gz" #coge la muestra pareada 
    #Vamos a comprobar que el archivo no está vacío
    { echo "Cheking if the "$sample" exists and it's not empty" 
    if [ -s "$file" ]; then
            echo "File exist and is not empty"  
    else
            echo "file is empty"  >&2
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastq
    echo "Checking whether file extension is .fastq.gz" 
    if [[ "$file" == *.fastq.gz ]]; then
            echo "File extension is correct" 
    else
            echo "File extension not correct. Must be .fastq.gz file" >&2
            exit 1
    fi 

    #Vamos a comprobar que el archivo es legible 
    echo "Checking "$sample" permissions" 

    if [[ -r "$file" ]]; then #solo necesito que sea legible
        echo "File is readable"
    else
        echo "File is not readable. Fixing." >&2
        chmod +r "$file"
    fi

    } 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out) 
    #Ahora vamos a trimar. 
    #Necesitamos saber como se llaman los archivos fastq y cómo se diferencian entre pares de secuencias
    #Va a sacar 4 archivos, dos de ellos pasarán la selección tras el trimado (paired) y otros dos no (unpaired). Esto ocurre porque son secuencias pareadas
    echo "Processing "$sample"..."
    #name=${sample%.*} #Para quitar la extensión y poder poner el nombre de la muestra en los archivos input y putput
    {
      java -jar trimmomatic-0.39.jar PE -phred33 \
      "$file1" "$file2" \
      "$output_dir/${sample}_1_paired.fastq.gz" "$output_dir/${sample}_1_unpaired.fastq.gz" \
      "$output_dir/${sample}_2_paired.fastq.gz" "$output_dir$/{sample}_2_unpaired.fastq.gz" \
      ILLUMINACLIP:$ILLUMINACLIP \
      SLIDINGWINDOW:$SLIDINGWINDOW \
      LEADING:$LEADING TRAILING:$TRAILING \
      MINLEN:$MINLEN
     } 2>> >(tee -a logs/${sample}_TRIM.err) >> >(tee -a logs/${sample}_TRIM.out) 
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
      fastqc "$output_dir/${sample}_1_paired.fastq.gz" "$output_dir/${sample}_2_paired.fastq.gz" -o "$fastqc_dir"
      } 2>> >(tee -a logs/fastqc/${sample}_fastqc.err)  >> >(tee -a logs/fastqc/${sample}_fastqc.out)
done


# Hacer multiqc para integrar el análisis del fastqc
{ echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$fastqc_dir" -n "multiqc_trimmed_analysis.html" -o "$multiqc_dir" } 2>> >(tee -a logs/multiqc/multiqc.err) >> >(tee -a logs/multiqc/multiqc.out)
   

echo -e "Analysis completed." >> >(tee -a logs/all_files.out)
###############################################
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
