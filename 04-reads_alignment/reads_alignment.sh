#!/bin/bash
# Este script se encarga de indexar el genoma de referencia y realizar los alineamientos.
# La herramienta empleada para el alineamiento es Hisat2:
# https://github.com/DaehwanKimLab/hisat2
# También se empleó samtools para convertir de .sam a .bam:
# https://www.htslib.org/
# Al final se realiza un MultiQC:
# https://github.com/MultiQC/MultiQC

#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.0"

# Usage example:
# ./alignment_script.sh -d input_reads_dir -G genome/genome.fa -A genome/annotation.gtf -o out_dir
#-v and -h are available for help and version

#######################################################################
# inicialización de logs vacíos
cat /dev/null > logs/
if [[ ! -e "logs/all_files.out" ]] && [[ ! -e "logs/all_files.err" ]]; then
    touch "logs/all_files.out"
    touch "logs/all_files.err"
fi

while getopts "hvd:G:A:o:" opt; do
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
        \?)  echo -e "Error: Invalid option. This script indexes the genome\n" \
           "Usage example:\n $0 -d input_reads_dir -G genome/genome.fa -A genome/annotation.gtf out_dir"
            echo "use -h for help and -v for version" >&2
            exit 1;;
    esac
done 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out) 

# Validación de archivos
# Convertir en bucle para dar permisos
for file in "$GENOME" "$GTF" "$INPUT_DIR"; do
    if [[ ! -e "$file" ]];then
        echo "Error: $file does not exist" >&2
        exit 1
    # Como los archivos se encuentran en una carpeta que no es nuestra, esto no tiene sentido porque no podemos hacer chmod
    # elif [ -r "$file" ]; then
    #     echo "File is readable but not executable. Fixing."
    #     chmod +x "$file"
    #     echo "fixed"
    # elif [ -x "$file" ]; then
    #     echo "File is executable but not readable"
    #     chmod +r "$file"
    #     echo "Fixed"
    # else
    #     echo "File is neither readable nor executable. Fixing"
    #     chmod +rx "$file"
    #     echo "Fixed"
    fi
done  2>> >(tee -a logs/${file}.err)  >> >(tee -a logs/${file}.out) 

#Creación del directorio:
echo "Creating output directory..."
if [[ ! -e "$OUTPUT_DIR" ]]; then
	echo "Output directory does not exists, creating..." | tee -a logs/output_dir.out
	mkdir "$OUTPUT_DIR"
fi 2>> >(tee -a logs/output_dir.err)  >> >(tee -a logs/output_dir.out) # los errores antes del 
# bucle importante van en otro archivo

# Indexación del genoma
echo "Running hisat2-build..." | tee -a logs/stdout
hisat2-build "$GENOME" "$OUTPUT_DIR/genome_index" 2>> >(tee -a all_files.out/err)

# Alineamiento
for file in "$INPUT_DIR"/*_1_filtered.fastq.gz; do
    sample=$(basename "$file" "_1_filtered.fastq.gz")
    file1="${INPUT_DIR}/${sample}_1_filtered.fastq.gz"
    file2="${INPUT_DIR}/${sample}_2_filtered.fastq.gz"
    
    # Crear carpeta de las muestras
    if [[ ! -e "$OUTPUT_DIR/$sample" ]]; then
    mkdir -p "$OUTPUT_DIR/$sample"
    fi
    
    echo "Aligning $sample ..." | tee -a logs/${sample}.out
    {
    # -x <genoma_ref> -1 <FW> -2 <RV> -S <output_sam_files> -p <threads>
    hisat2 -x "$OUTPUT_DIR/genome_index" \
    -1 "$file1" \
    -2 "$file2" \
    -S "$OUTPUT_DIR/$sample/HISAT2.sam" \
    --summary-file "$OUTPUT_DIR/$sample/${sample}_hisat2_summary.txt" \
    -p 10 && echo "Alignment with sample $sample done" || echo "Alignment with sample $sample failed"
    } 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out) 
    
    echo "Converting $sample to a sorted bam file..." | tee -a logs/${sample}.out
    {
    samtools view -bS "$OUTPUT_DIR/$sample/HISAT2.sam" > "$OUTPUT_DIR/$sample/${sample}.bam" && echo "Converting sample $sample done" || echo "Converting sample $sample failed"
    samtools sort "$OUTPUT_DIR/$sample/${sample}.bam" -o "$OUTPUT_DIR/$sample/${sample}_sorted.bam" && echo "Sorting sample $sample done" || echo "Sorting sample $sample failed"
    # samtools permite conviertir .sam en .bam, ocupa menos al ser los binarios y son necesarios para después
    # No nos interesa el sam, pero sí el sorted_bam para el conteo y el bam sirve también para hacer después un MultiQC.
    rm "$OUTPUT_DIR/$sample/HISAT2.sam"
    # rm "$OUTPUT_DIR/$sample/HISAT2.bam"
    } 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out)
    echo "Finished with $sample"

done

# añadir un MultiQC en una carpeta llamada /results/MultiQC
# Creación del MultiQC results
if [[ ! -e "$OUTPUT_DIR/multiqc" ]]; then
    mkdir -p "$OUTPUT_DIR/multiqc"
fi

# MultiQC con los resultados del alineamient
{
echo "Running MultiQC..." | tee -a logs/stdout
multiqc "$OUTPUT_DIR" -n "multiqc_alignment_analysis.html" -o "$OUTPUT_DIR/multiqc" --force
} 2>> >(tee -a logs/multiqc/multiqc.err) >> >(tee -a logs/multiqc/multiqc.out)
