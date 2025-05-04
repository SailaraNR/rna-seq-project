
#Comprobar que el archivo es .fastqc

while getopts "hv" opt; do
        case $opt in
                h) echo "You asked for usage help"
                   echo "This script will trim the raw seqs" 
                   echo -e "You won't need to add any flags, but -v and -h are available to ask for script version \n or help" ;;
                v) echo "Version 1.0" 
                   exit 1 ;;
                ?) echo "Invalid option or missing argument" >&2
                   echo "Usage example: $0 -u URL" >&2
                   exit 1;;
        esac
done

for files in ../01-pre_fastqc; do
    #Vamos a comprobar que el archivo es "readable" y ejecutable
    echo "Cheking $file permissions"
    if [ -r $file ]; then
        if [ -x $file ]; then
                echo "This file is readable and executable"
        else
                echo "This file is readable but not executable"
                chmod +x $file
                echo "Now is ready"
    elif [ -x $file ]; then
        echo "This file is not readable but executable"
        chmod +r $file
        echo "Now is ready"
    else
        echo "This file es not readable not executable"
        chmod +rx $file
        echo "Now is ready"

    fi
    #También se comprobará si tiene la extensión correcta (.fastqc)
    echo "Checking whether file extension is .fastq"
    if [ "$file" == *.fastq ]; then
        echo "File extension is correct"
    else
        echo "File extension not correct. Must be .fastq file" >&2
        exit 1
    fi
    #También se comprobará si el archivo está vacío
    echo "Cheking if the file exists and it's not empty"
    if [ -s $file ]; then
        echo "File exist and is not empty" 
    else
        echo "$file is empty"  >&2
        exit 1
    fi 
    #Ahora vamos a trimar

    trimmomatic PE \
      input_R1.fastq input_R2.fastq \
      output_R1_paired.fastq output_R1_unpaired.fastq \
      output_R2_paired.fastq output_R2_unpaired.fastq \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:4:20 \
      LEADING:3 TRAILING:3 \
      MINLEN:36
done