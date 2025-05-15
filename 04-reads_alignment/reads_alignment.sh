#!/bin/bash
# Este script se encarga de indexar el genoma de referencia para los alineamientos.
# La herramienta empleada para el alineamiento es Hisat2:
# https://github.com/DaehwanKimLab/hisat2
# También se empleó samtools para convertir de .sam a .bam:
# https://www.htslib.org/

#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.0"

# Usage example:
# ./alignment_script.sh -d input_reads_dir -G genome/genome.fa -A genome/annotation.gtf -o out_dir
#-v and -h are available for help and version

#######################################################################
# inicialización de logs vacíos
cat /dev/null > logs/*-v

while getopts "hvd:G:A:o:"; do
    case $opt in
        h) echo -e "This script indexes the genome\n"\
           "Usage example:\n $0 -d input_reads_dir -G genome/genome.fa -A genome/annotation.gtf -o out_dir" | tee -a logs/stderr
            echo "use -h for help and -v for version" | tee -a logs/stderr
            exit 1;;
        v) echo "$version" | tee -a logs/stdout
            exit 0;;
        d) INPUT_DIR="$OPTARG";;
        G) GENOME="$OPTARG";;
        A) GTF="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        \?)  echo -e "Error: Invalid option. This script indexes the genome\n"\
           "Usage example:\n $0 -d input_reads_dir -G genome/genome.fa -A genome/annotation.gtf out_dir"
            echo "use -h for help and -v for version" >&2
            exit 1;;
    esac
done

# Validación de archivos
# Convertir en bucle para dar permisos
for file in "$GENOME" "$GTF" "$INPUT_DIR"; do
    if [[ ! -e "$file" ]];then
        echo "Error: $file does not exist" | tee -a logs/${file}.err
        exit 1
    elif [[ -r "$file" && -x "$file" ]]; then
        echo "File is readable and executable" | tee -a logs/${file}.out
    elif [ -r "$file" ]; then
        echo "File is readable but not executable. Fixing." | tee -a logs/${file}.out
        chmod +x "$file"
        echo "fixed" | tee -a logs/${file}.out
    elif [ -x "$file" ]; then
        echo "File is executable but not readable" | tee -a logs/${file}.out
        chmod +r "$file"
        echo "Fieexd" | tee -a logs/${file}.out
    else
        echo "File is neither readable nor executable. Fixing" | tee -a logs/${file}.out
        chmod +rx "$file"
        echo "Fixed" | tee -a logs/${file}.out
    fi
done  2>> >(tee -a logs/${file}.err)  >> >(tee -a logs/${file}.out) 

#Creación del directorio:
echo "Creating output directory..."
if ! [[ -e "$OUTPUT_DIR" ]]; then
	echo "Output directory does not exists, creating..." | tee -a logs/output_dir.out
	mkdir "$OUTPUT_DIR"
fi


# Indexación del genoma
echo "Running hisat2-build..." | tee -a logs/stdout
hisat2-build "$GENOME" "$OUTPUT_DIR/genome_index" 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out) 

# Alineamiento
for file in "$INPUT_DIR"/*_1_filtered.fastq.gz; do
    sample=$(basename "$file" _1_filtered.fastq.gz)
    file1="${INPUT_DIR}/${sample}_1_filtered.fastq.gz"
    file2="${INPUT_DIR}/${sample}_2_filtered.fastq.gz"
    
    # Crear carpeta de las muestras
    if [[ ! -e "$OUTPUT_DIR/$sample/results"]]; then
    mkdir -p "$OUTPUT_DIR/$sample/results"
    fi
    
    echo "Aligning $sample ..." | tee -a logs/${sample}.out
    {
    # -x <genoma_ref> -1 <FW> -2 <RV> -S <output_sam_files> -p <threads>
    hisat2 -x "$OUTPUT_DIR/genome_index" \
    -1 "$file1" \
    -2 "$file2" \
    -S "$OUTPUT_DIR/$sample/results/HISAT2.sam" \
    -p 6 && echo "Alignment with sample $sample done" || echo "Alignment with sample $sample failed"

    samtools view -bS "$OUTPUT_DIR/$sample/results/HISAT2.sam" > "$OUTPUT_DIR/$sample/results/HISAT2.bam"
    samtools sort "$OUTPUT_DIR/$sample/results/HISAT2.bam" -o "$OUTPUT_DIR/$sample/results/${sample}_sorted.bam"
    # samtools permite conviertir .sam en .bam, ocupa menos al ser los binarios y son necesarios para después
    # Se borran el .sam y el .bam desordenado. Lo único que nos interesa es el sorted .bam
    rm "$OUTPUT_DIR/$sample/results/HISAT2.sam"
    rm "$OUTPUT_DIR/$sample/results/HISAT2.bam"

    } 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out) 
    echo "Finished $sample" | tee -a logs/${sample}.out
done

# añadir un MultiQC en una carpeta llamada /results/MultiQC
# Creación del MultiQC results
if [[ ! -d "$OUTPUT_DIR/MultiQC" ]]; then
    mkdir -p "$OUTPUT_DIR/MultiQC"
fi

# MultiQC con los resultados del alineamient
{
echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$OUTPUT_DIR" -n "multiqc_alignment_analysis.html" -o "$OUTPUT_DIR/MultiQC" \
} 2>> >(tee -a logs/multiqc/multiqc.err) >> >(tee -a logs/multiqc/multiqc.out)