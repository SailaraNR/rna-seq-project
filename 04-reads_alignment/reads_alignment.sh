# !/bin/bash
# Este script se encarga de indexar el genoma de referencia para los alineamientos.
# La herramienta empleada es Hisat2 cuyo GitHub es el siguiente:
# https://github.com/DaehwanKimLab/hisat2
#Author: Laura Barreales y Sara Lévano
#Start date: 10st may 2025
version="versión 1.0"

# Usage example:
# ./alignment_script.sh -d input_reads_dir -F _R1.fastq -R _R2.fastq -G genome/genome.fa -A genome/annotation.gtf -t STAR out_dir -l sample_list.txt
#-v and -h are available for help and version

#######################################################################

while getops "hv:d:F:R:G:A:t:o:l"; do
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
if [[ ! -e "$GENOME" ]] || [[ ! -r "$GENOME" ]];then
	echo "Error: $GENOME does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
elif [[ ! -e "$GTF" ]] || [[ ! -r "$GTF" ]]; then
	echo "Error: $GTF does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
elif [[ ! -e "$sample_list" ]] || [[ ! -r "$sample_list" ]]; then
	echo "Error: $sample_list does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
elif [[ ! -e "$INPUT_DIR" ]] || [[ ! -r "$INPUT_DIR" ]]; then
	echo "Error: $INPUT_DIR does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
fi

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
    {
    # Vale, aquí tengo dudas porque se puede usar también con --sra-acc <SRA accession number>, -q para .fasta, etc. y no sé qué es más práctico
    # -x <genoma_ref> -1 <FW> -2 <RV> -S <output_sam_files>
    hisat2 -x $OUTPUT_DIR/HISAT2_genome_indices/genome_index \
    -1 $INPUT_DIR/$FW \
    -2 $INPUT_DIR/$RV \
    -S $OUTPUT_DIR/$sample/results/HISAT2/HISAT2.sam && echo "Alignment with sample $sample done" || echo "Alignment with sample $sample failed"
    -p/--threads 3
    } >> logs/stdout 2>> logs/stderr
done < sample_list