#!/bin/bash
#Este script coge los archivos con las raw data, las analiza, filtra y elimina secuencias de baja calidad
#Los archivos que acepta solo son archivos .fastq.gz 
#Author: Laura Barreales and Sara Lévano
#Start date: 14th May 2025
version="Version: 1.3"

#Usage example: ./fastp.sh


# Directorios que deberían ponerse
#input_dir="../00-raw_data/results"
#output_dir="./results"
#fastqc_dir="./results/fastqc"
#multiqc_dir="./results/multiqc"

#TENEMOS QUE CREAR UNA CARPETA PARA MULTIQC Y FASTQC DENTRO DE LOS LOGS
#######################################################################
# inicialización de logs vacíos
cat /dev/null > logs/*

#Usage
{ echo -e "Use -h for help or -v for the version. \nIntroudce input, output, fastqc and multiqc directories, but do not introduce more arguments unless you want to change filter parameters." 

    # Parámetros por defecto para fastp por si el usuario no pone ninguno
    #Si el usuario pone alguno, estas variables se reescribiran gracias al getopts 
    MIN_QUAL=20
    MIN_LEN=40 #ente 36-50 suele estar bien

    #Manejar opciones y argumentos
    while getopts "hvm:l:i:o:f:c:" opt; do
        case $opt in
            h) echo -e "Usage help:\nThis script filters and trims raw sequences using fastp."
            echo -e "Usage:\n\t$0 [-m MIN_QUAL] [-l MIN_LEN] [-i INPUT] [-o OUTPUT] [-f FASTQC_DIR] [-c MULTIQC_DIR]"
            exit 0;;
            v) echo "$version"
            exit 0;;
            m) MIN_QUAL="$OPTARG" ;;
            l) MIN_LEN="$OPTARG" ;;
            i) input_dir="$OPTARG" ;;
            o) output_dir="$OPTARG" ;;
            f) fastqc_dir="$OPTARG" ;;
            c) multiqc_dir="$OPTARG" ;;
            \?) echo -e "Invalid option\nUsage: $0 ..." >&2
                exit 1;; #si pone simbolos raros no sirve
        esac
    done
    echo "Using minimun quality = $MIN_QUAL and minimun length = $MIN_LEN"



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
} 2>> >(tee -a logs/all_files.err) >> >(tee -a logs/all_files.out)

# Procesar todos los archivos .fastq.gz
for file in "$input_dir"/*_1.fastq.gz; do #Los archivos fastq para trimar y filtrar se encuentran en esa carpeta. $file debería ser una ruta ansoluta al file
    sample=$(basename "$file" "_1.fastq.gz") #Nos quedamos con el nombre de la muestra. No extensión
    file1="${input_dir}/${sample}_1.fastq.gz"
    file2="${input_dir}/${sample}_2.fastq.gz" #coge las muestras pareadas
    #Vamos a comprobar que el archivo no está vacío
    { echo "Checking if $sample exists and is not empty"
    if [ -s "$file1" ] && [ -s "$file2" ]; then
        echo "Files exist and are not empty"
    else
        echo "One or both paired files are empty or missing" >&2
        exit 1
    fi
    #Vamos a comprobar que el arhivo es .fastq.gz
    echo "Checking whether file extension is .fastq.gz"
    if [[ "$file1" == *.fastq.gz && "$file2" == *.fastq.gz ]]; then 
        echo "File extensions are correct"
    else
        echo "Incorrect file extensions" >&2
        exit 1
    fi
    #Vamos a comprobar que ambos archivos son legibles
    echo "Checking read permissions"
    for f in "$file1" "$file2"; do
        if [[ ! -r "$f" ]]; then #solo necesito que sean legibles
            echo "File $f not readable. Attempting to fix." >&2
            chmod +r "$f"
        fi
    done
    } 2>> >(tee -a logs/${sample}.err) >> >(tee -a logs/${sample}.out)
    #Ahora vamos a trimar
    echo "Processing $sample..."

    {
    fastp \
        -i "$file1" -I "$file2" \
        -o "$output_dir/${sample}_1_filtered.fastq.gz" \
        -O "$output_dir/${sample}_2_filtered.fastq.gz" \
        -q "$MIN_QUAL" -l "$MIN_LEN" \
        --detect_adapter_for_pe \
        --thread 4 \
        --html "$output_dir/${sample}_fastp.html" \
        --json "$output_dir/${sample}_fastp.json" \
        --trim_poly_x --poly_x_min_len 10 \
        --trim_front1 10 --trim_front2 10\
        --report_title "$sample fastp report"  
        #i: input files, son los archivos del 00/results
        #o: outputfiles, serán los archivos de resultados que se guardan en el outputdir
        #q: calidad mínima para trimar o no
        #l: logitud minima después del trimado
        #detecta adaptadores
        #usa 4 hilos
        #genera informes  .html y .json
        #detecta polyA y polyT que sean de mínimo 10 bases y los elimina
        #el titulo del report 
    } 2>> >(tee -a logs/${sample}_fastp.err) >> >(tee -a logs/${sample}_fastp.out)
#ahora vamos a ejecutar fastqc para las seqs trimadas, están en 01/results y los guarda en 01/results/fastqc
    {
    echo "Executing FastQC for $sample" 
    fastqc "$output_dir/${sample}_1_filtered.fastq.gz" "$output_dir/${sample}_2_filtered.fastq.gz" -o "$fastqc_dir"
    } 2>> >(tee -a logs/fastqc/${sample}_fastqc.err) >> >(tee -a logs/fastqc/${sample}_fastqc.out)
done
#ahroa ejecutamos multiqc cogiendo los resultados del fasqc del 01/reults/fastqc y los guarda en 01/results/multiqc
{ echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$fastqc_dir" -n "multiqc_fastp_analysis.html" -o "$multiqc_dir"
} 2>> >(tee -a logs/multiqc/multiqc.err) >> >(tee -a logs/multiqc/multiqc.out)

echo -e "Analysis completed." >> >(tee -a logs/all_files.out)
