# !/bin/bash
# Este script se encarga de contar el múmero de lecturas de un alineamiento partiendo de archivos .bam
# La herramienta usada es FeatureCounts:
# https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html

#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.0"

# Usage example:
# ./counts_alignment.sh -A genome_annotation.gtf -o output.txt -d dir_sorted_example_alignment
#-v and -h are available for help and version

#######################################################################
#inicialización de logs
cat /dev/null > logs/
if [[ ! -e "logs/all_files.out" ]] && [[ ! -e "logs/all_files.err" ]]; then
    mkdir -p  "logs/all_files.out"
    mkdir -p  "logs/all_files.err"
fi

while getops "hvA:o:f:"; do
    case $opt in
        v) echo "$version" | tee -a logs/stdout
        exit 0;;
        h) echo -e "This script counts the number of reads of an aligment using a bam file and  GTF\nUsage example:\n"\
        "./counts_alignment.sh -A genome_annotation.gtf -o output.txt -f sorted_example_alignment.bam" | tee -a logs/stdout
        exit 0;;
        A) GTF="$OPTARG";;
        o) output="$OPTARG";;
        d) alignment="$OPTARG";;
        \?) echo -e "Invalid option: This script counts the number of reads of an aligment using a bam file and  GTF\nUsage example:\n"\
        "./counts_alignment.sh -A genome_annotation.gtf -o output.txt -d dir_sorted_example_alignment" >&2
        exit 1;;
    esac
done 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out) 

# File validation: Check if the variables exist
# El GTF no es nuestro, así que tampoco podemos cambiar los permisos, solo verificar si existe o no
if [[ ! -e "$GTF" ]] || [[ ! -r "$GTF" ]]; then
	echo "Error: $GTF does not exist or does not have read permissions" >&2
	exit 1
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out) 
if [[ ! -e "$alignment" ]] || [[ ! -r "$alignment" ]]; then
	echo "Error: $alignment does not exist or does not have read permissions" >&2
	exit 1
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out) 

#Create output file
echo "Creating output.txt file..."
if ! [[ -e $output ]]; then
	echo "Output file does not exists, creating..."
	touch $output
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)
{
for file in "$alignment"/SRR*/*.bam; do # 04-reads_alignment/results en nuestro caso sería el $alignment
        sample=$(basename "$file" "_sorted.bam")
        echo -e "\nChecking if $sample exists and it's not empty"
        if [[ -f "$file" && -s "$file" ]]; then
            echo "$sample exists and is not empty"
        else
            echo "$sample is not a file or is empty"  >&2
            exit 1
        fi
        echo "Checking $sample permissions"
            if [[ ! -r "$file" ]]; then
                echo "File is not readable. Fixing..." >&2
                chmod +r "$file" && echo "Permissions fixed" || echo "Could not change permisions" >&2
                    continue
            fi
    done

    #Counting reads
    featureCounts -p -O -T 6 -a "$GTF" -o "$output" "$file" && echo "Counts for "$sample" done" || echo "Counts for "$sample" failed"
} 2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out)

# -p species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
# -t type of fragments (exon, intron...)
# -O assigns reads to all their overlapping meta-features.
# -T specifies the number (6) of threads to be used.
# -a is the genome annotation file ($GTF).
# -o specifies the name of the output file, which includes the read counts ($output).
# $alignment is an alignment file: in this file, the reads we want to count are aligned to the same genome as the annotation file