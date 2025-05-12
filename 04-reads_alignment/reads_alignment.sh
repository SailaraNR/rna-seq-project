# !/bin/bash
# Este script se encarga de indexar el genoma de referencia para los alineamientos.
# La herramienta empleada para el alineamiento es Hisat2:
# https://github.com/DaehwanKimLab/hisat2
# También se empleó samtools para convertir de .sam a .bam:
# https://www.htslib.org/

#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.0"

# Usage example:
# ./alignment_script.sh -d input_reads_dir -F _R1.fastq -R _R2.fastq -G genome/genome.fa -A genome/annotation.gtf -t STAR out_dir -l sample_list.txt
#-v and -h are available for help and version

#######################################################################

while getops "hvd:F:R:G:A:t:o:l:"; do
    case $opt in
        h) echo -e "This script indexes the genome\n"\
           "Usage example:\n $0 -d input_reads_dir -F _R1.fastq -R _R2.fastq -G genome/genome.fa -A genome/annotation.gtf -t STAR out_dir -l sample_list.txt" | tee -a logs/stderr
            echo "use -h for help and -v for version" | tee -a logs/stderr
            exit 1;;
        v) echo "$version" | tee -a logs/stdout
            exit 0;;
        d) INPUT_DIR="$OPTARG";;
        F) FW="$OPTARG";;
        R) RV="$OPTARG";;
        G) GENOME="$OPTARG";;
        A) GTF="$OPTARG";;
        t) align_tool="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        l) sample_list="$OPTARG";;
        \?)  echo -e "Error: Invalid option. This script indexes the genome\n"\
           "Usage example:\n $0 -d input_reads_dir -F _R1.fastq -R _R2.fastq -G genome/genome.fa -A genome/annotation.gtf -t STAR out_dir -l sample_list.txt" | tee -a logs/stderr
            echo "use -h for help and -v for version" | tee -a logs/stderr
            exit 1;;
    esac
done

# File validation: Check if the variables exist
# Convertir en bucle para dar permisos
for file in "$GENOME" "$GTF" "$sample_list" "$INPUT_DIR"; do
    if [[ ! -e "$file" ]];then
        echo "Error: $file does not exist" | tee -a logs/stderr >&2
        exit 1
    elif [[ -r "$file" && -x "$file" ]]; then
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
done >> logs/stdout 2>> logs/stderr

#Create output directory and subdirectories:
echo "Creating output directory..."
if ! [[ -e $OUTPUT_DIR ]]; then
	echo "Output directory does not exists, creating..." | tee -a logs/stout
	mkdir $OUTPUT_DIR
fi

#Index creation:
hisat2-build $genome_fasta $OUTPUT_DIR/HISAT2_genome_indices/genome_index

# Aligment
while IFS= read -r sample; do
    echo "Aligning $sample ..." | tee -a logs/stout
    {
    # -x <genoma_ref> -1 <FW> -2 <RV> -S <output_sam_files>
    hisat2 -x $OUTPUT_DIR/HISAT2_genome_indices/genome_index \
    -1 $INPUT_DIR/$FW \
    -2 $INPUT_DIR/$RV \
    -S $OUTPUT_DIR/$sample/results/HISAT2/HISAT2.sam && echo "Alignment with sample $sample done" || echo "Alignment with sample $sample failed"
    samtools view -bS $OUTPUT_DIR/$sample/results/HISAT2/HISAT2.sam > $OUTPUT_DIR/$sample/results/HISAT2/HISAT2.bam 
    samtools sort $OUTPUT_DIR/$sample/results/HISAT2/HISAT2.bam -o $OUTPUT_DIR/$sample/results/HISAT2/sorted_HISAT2.bam
    # samtools permite conviertir .sam en .bam, ocupa menos al ser los binarios
    # -p/--threads 3 # En caso de querer especificar los hilos
    } >> logs/stdout 2>> logs/stderr
    echo "Finished $sample" | tee -a logs/stout
done < sample_list

# añadir un MultiQC en una carpeta llamada /results/MultiQC
# command 2>&1 >file (mirar esta wea de las redirecciones)