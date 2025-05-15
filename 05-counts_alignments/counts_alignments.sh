# !/bin/bash
# Este script se encarga de contar el múmero de lecturas de un alineamiento partiendo de archivos .bam
# La herramienta usada es FeatureCounts:
# https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html

#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.0"

# Usage example:
# ./counts_alignment.sh -A genome_annotation.gtf -o output.txt -f sorted_example_alignment.bam
#-v and -h are available for help and version

#######################################################################

while getops "hvA:o:f:"; do
    case $opt in
        v) echo "$version" | tee -a logs/stdout
        exit 0;;
        h) echo -e "This script counts the number of reads of an aligment using a bam file and  GTF\nUsage example:\n"\
        "./counts_alignment.sh -A genome_annotation.gtf -o output.txt -f sorted_example_alignment.bam" | tee -a logs/stdout
        exit 0;;
        A) GTF="$OPTARG";;
        o) output="$OPTARG";;
        f) alignment="$OPTARG";;
        \?) echo -e "Invalid option: This script counts the number of reads of an aligment using a bam file and  GTF\nUsage example:\n"\
        "./counts_alignment.sh -A genome_annotation.gtf -o output.txt -f sorted_example_alignment.bam" | tee -a logs/stderr
        exit 1;;
    esac
done

# File validation: Check if the variables exist
if [[ ! -e "$GTF" ]] || [[ ! -r "$GTF" ]]; then
	echo "Error: $GTF does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
fi
if [[ ! -e "$alignment" ]] || [[ ! -r "$alignment" ]]; then
	echo "Error: $alignment does not exist or does not have read permissions" | tee -a logs/stderr
	exit 1
fi

#Create output file
echo "Creating output.txt file..."
if ! [[ -e $output ]]; then
	echo "Output file does not exists, creating..." | tee -a logs/stout
	touch $output
fi

#Counting reads
featureCounts -p -t exon -O -T 6 -a "$GTF" -o "$output" "$alignment" | tee -a logs/stdout 2>> logs/stderr

# -p species that fragments (or templates) will be counted instead of reads. This is only applicable for paired-end reads.
# -t type of fragments (exon)
# -O assigns reads to all their overlapping meta-features.
# -T specifies the number (6) of threads to be used.
# -a is the genome annotation file ($GTF).
# -o specifies the name of the output file, which includes the read counts ($output).
# $alignment is an alignment file: in this file, the reads we want to count are aligned to the same genome as the annotation file