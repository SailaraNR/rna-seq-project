#Este script coge los archivos con las raw data, las analiza, trima y elimina secuencias
#Los archivos que acepta solo son archivos .fastq 
#Author: Laura Barreales and Sara Lévano
#Start date: 6th May 2025
#Version: 1.0

#Usage example: trimmingfiltering.sh
#No arguments expected, but -v and -h are available

#######################################################################

#explain the code
echo -e "Use -h for help or -v for the version. \nBut do not introduce any arguments or options for for trimming and filtering sequences" | tee -a logs/stdout

#Manejar opciones y argumentos
while getopts "hv" opt; do
    case $opt in
        h) echo -e "You asked for usage help\n"\
            "This script trims and filters raw sequences" | tee -a logs/stdout
            echo -e "Usage example:\n\t sh $0 -h: Provides help \n\t sh $0 -v: Tells script's version\n\tsh $0: trims and filters raw sequences" | tee -a logs/stdout
            exit 0;;
        v) echo "Version 1.0" | tee -a logs/stdout
            exit 0;;
        \?) echo -e "Invalid option\n\ 
            Usage example: $0 " | tee -a logs/stdout
            exit 1;; #si pone simbolos raros no sirve
    esac
done

for file in ../01-pre_fastqc/results; do
    #Vamos a comprobar que el archivo no está vacío
    echo "Cheking if the file exists and it's not empty" | tee -a logs/stdout
    if [ -s $file ]; then
            echo "File exist and is not empty"  | tee -a logs/stdout
    else
            echo "file is empty"  | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el arhivo es .fastq
    echo "Checling whether file extension is .fastq" | tee -a logs/stdout
    if [ "$file" == *.fastq ]; then
            echo "File extension is correct" | tee -a logs/stdout
    else
            echo "File extension not correct. Must be .txt file" | tee -a logs/stderr
            exit 1
    fi 

    #Vamos a comprobar que el archivo es legible y ejecutable
    echo "Cheking $file permissions" | tee -a logs/stdout

    if [[ -r "$file" && -x "$file" ]]; then
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

    #Ahora vamos a trimar. 
    #Necesitamos saber como se llaman los archivos fastq y cómo se diferencian entre pares de secuencias
    #Va a sacar 4 archivos, dos de ellos pasarán la selección tras el trimado (paired) y otros dos no (unpaired). Esto ocurre porque son secuencias pareadas
    echo "Processing $file..."
    trimmomatic PE \
      {$file}_R1.fastq {$file}_R2.fastq \
      output_R1_paired.fastq output_R1_unpaired.fastq \
      output_R2_paired.fastq output_R2_unpaired.fastq \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:4:20 \
      LEADING:3 TRAILING:3 \
      MINLEN:36 
      #MINLEN: Es la cantidad mínima de bases que pedimos tras el trimado, yo pondría más
      #LEADING: Es la calidad mínima que debe tener la secuenciación de cada base del principio, si no llega se corta
      #TRAILING: Es la calidad mínima que debe tener la secuenciación de cada base del final, si no llega se corta
      #SLIDINGWINDOW:<windowSize>:<requiredQuality>  Mira en la ventana de bases cual es la calidad media, si no supera un mínimo corta. Compatible con TRAILING y LEADING?
      #ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>  Preguntar si esas secuencias tienen adaptadores
      #TOPHRED33/64 This (re)encodes the quality part of the FASTQ file to base 33 or 64.
      #NO utilizaría MAXINFO  
done