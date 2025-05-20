# !/bin/bash
# Este script se encarga de contar el múmero de lecturas de un alineamiento par$
# La herramienta usada es FeatureCounts:
# https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html
#Author: Laura Barreales y Sara Lévano
#Start date: 10th may 2025
version="versión 1.3"
# Usage example:
# ./counts_alignment.sh -A genome_annotation.gtf -o output.txt -d dir_sorted_ex$
#-v and -h are available for help and version
#######################################################################

while getopts "hvA:o:d:" opt; do
    case $opt in
        v) echo "$version" 
        exit 0;;
        h) echo -e "This script counts the number of reads of an aligment using ./$0 -A genome_annotation.gtf -o output.txt -f sorted"
        exit 0;;
        A) GTF="$OPTARG";;
        o) output="$OPTARG";;
        d) alignment="$OPTARG";;
        \?) echo -e "Invalid option: This script counts the number of reads of \n./$0 -A genome_annotation.gtf -o output.txt -d dir_so"
        exit 1;;
    esac
done 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

# File validation: Check if the variables exist
# El GTF no es nuestro, así que tampoco podemos cambiar los permisos, solo verificar si existe o no
if [[ ! -e "$GTF" ]] || [[ ! -r "$GTF" ]]; then
        echo "Error: $GTF does not exist or does not have read permissions" >&2
        exit 1
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

#Comprobar si existe 04/results (en principio debería existir porque se lo damos predetermiando
#pero podemos mirar si está vacío)
{
    if [[ ! -e "$alignment" ]]; then
        echo "Error: $alignment does not exist" >&2
        exit 1
    fi

    if [[ ! -r "$alignment" ]]; then
        echo "$alignment is not readable. Fixing"
        chmod +r "$alignment"
        echo "Fixed"
    fi

    if [[ ! -x "$alignment" ]]; then
        echo "$alignment is not executable. Fixing"
        chmod +x "$alignment"
        echo "Fixed"
    fi
} 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)

#Create output file
echo "Creating output.txt file..."
if ! [[ -e $output ]]; then
        echo "Output file does not exists, creating..."
        touch $output
fi 2>> >(tee -a logs/all_files.err)  >> >(tee -a logs/all_files.out)


for file in ${alignment}/SRR*/*.bam; do # 04-reads_alignment/results en nuestro caso sería el $alignment
    sample=$(basename "$file" "_sorted.bam")
    {
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
    #Counting reads
    featureCounts -p -O -T 6 -a "$GTF" -o "$output" "$file" && echo "Counts for "$sample" done" || echo "Coun"
    }  2>> >(tee -a logs/${sample}.err)  >> >(tee -a logs/${sample}.out)
done 

echo "Counting has finished"
# -p species that fragments (or templates) will be counted instead of reads. This is only applicable for pair$
# -t type of fragments (exon, intron...)
# -O assigns reads to all their overlapping meta-features.
# -T specifies the number (6) of threads to be used.
# -a is the genome annotation file ($GTF).
# -o specifies the name of the output file, which includes the read counts ($output).
# $file in $alignment is an alignment file: in this file, the reads we want to count are aligned to the same 

