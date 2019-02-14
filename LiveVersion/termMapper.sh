#!/usr/bin/env bash
#Package Version: 4.1

############################################################################################################
# Author(s): T.J.Cooper
# Updated: 15/4/2018
# Automates batch-processing of FASTQ/SAM files for genome-wide mapping of terminally informative molecules
############################################################################################################

function usage {
   echo -e "\nUsage: termMapper -i [INPUT FOLDER] -c [CONFIGURATION FILE] -o [OUTPUT FOLDER]\n";
   echo -e "-i INPUT: Input data folder containing paired-end FASTQ files.\n"
   echo -e "-c CONFIG: Configuration file specifying user-parameters (termMapper.config).\n";
   echo -e "-o OUTPUT: Output data folder.\n"; exit 1;
}

###########################
# Dependencies
###########################
command -v perl >/dev/null 2>&1 || { echo -e >&2 "\nError: Perl installation not found.\n"; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo -e >&2 "\nError: Bowtie2 installation not found.\n"; exit 1; }
[ "${BASH_VERSINFO:-0}" -lt 4 ] && { echo -e "\nError: Bash (v4.0+) installation not found.\n"; exit 1; }
BOWTIE="$(bowtie2 --version)"
BTVER=$(echo "$BOWTIE" | awk '/version/ {print $NF}')

###########################
# Command-line Options
###########################
CONF=""
INPUT=""
RES=""
while getopts ":c:i:o:" FLAG; do
   case $FLAG in
      i) INPUT=$OPTARG;;
      c) CONF=$OPTARG;;
      o) RES=$OPTARG;;
      \?) echo -e "\nInvalid option: -$OPTARG"
         usage;;
      :) echo -e "\nOption -$OPTARG requires an argument."
         usage;;
   esac
done
if [ "$#" -eq 0 ]; then
    usage
fi

#################################
# Read Config File & Directories
#################################
BASEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
PERLDIR=$BASEDIR/Scripts/readProcessing.pl
shopt -s extglob
declare -A config=()
while IFS='=' read -r k v; do
   [[ $v ]] || continue
   [[ $k = "#"* ]] && continue
   k=${k%+([[:space:]])}
   v=${v#+([[:space:]])}
   config[$k]=$v
done < "$CONF"
if [[ ! -f $PERLDIR || -z ${config[READ1_EXT]} || -z ${config[READ2_EXT]} || -z ${config[CORE]} || -z ${config[SPACE_SAVER]} ]]; then
    echo -e "\nError (Config File): User parameters or script files are missing and/or incorrectly specified.\n"; exit -1
fi
if [[ -z ${config[MAPQ_FILTER]} ]]; then
    echo -e "\nError (Config File): MAPQ threshold unspecified (Off = 0, Max = 42).\n"; exit -1
fi
if [[ ! -d ${config[GENOME_DIR]} || -z ${config[GENOME_NAME]} || ! -f ${config[FASTA_INDEX]} ]]; then
    echo -e "\nError (Config File): Genome files/directories are incorrectly specified.\n"; exit -1
fi
if [[ -d "$RES" && "$(ls -A $RES)" ]]; then
  while true; do
    read -p "$(echo -e '\nSpecified results folder already contains files. Do you wish to overwrite?\n\b')" yn
    case $yn in
      [Yy]* ) rm -rf $RES/*; break;;
      [Nn]* ) exit;;
      * ) echo "Please answer Yes or No.";;
    esac
  done
fi
mkdir -p "$RES" || exit
mkdir -p "$RES/Logs" || exit
exec > >(tee "$RES/Logs/System.log")

###########################
# Input Files
###########################
cd "$INPUT" || exit
declare -A STR
shopt -s nullglob
for file in *${config[READ1_EXT]}.*; do
   STR["${file%*${config[READ1_EXT]}.*}"]=1
done
if [ ${#STR[@]} = 0 ]; then
   echo -e "\nNo valid FASTQ files were found. Please verify that the READ1_EXT/READ2_EXT variables are correctly set (termMapper.config). "
   usage
fi
cp $CONF $RES/Logs/

###########################
# Bowtie2 Alignment
###########################
cd "$RES" || exit
ST=$(date +%s)
echo "-------------------------------------"
echo "FASTQ Alignment (--end-to-end)"
echo "-------------------------------------"
echo "Currently aligning:"
for i in "${!STR[@]}"; do
  printf "%s\n%s\n%s\n%s\n" "-----------------------" "Bowtie2 Version: $BTVER" "Reference Genome: ${config[GENOME_NAME]}" "Sample: $i" > "$RES/Logs/Log.$i.txt"
   echo ${i}
   printf "%s\n%s\n%s\n%s\n" "---------------------------------------------------------------------" "GLOBAL" "Parameters: ${config[GLOBAL_OPTIONS]}" "---------------------------------------------------------------------" >> "$RES/Logs/Log.$i.txt"
   bowtie2 --end-to-end ${config[GLOBAL_OPTIONS]} -p ${config[CORE]} -x ${config[GENOME_DIR]}/${config[GENOME_NAME]} -1 $INPUT/$i${config[READ1_EXT]}.fastq -2 $INPUT/$i${config[READ2_EXT]}.fastq -S $RES/$i"_Global".SAM 2>> "$RES/Logs/Log.$i.txt"
done
AT=$(date +%s)









###########################
# Data Processing
###########################
echo "-------------------------------------"
echo "Calculating Coordinates....";
echo "-------------------------------------"
echo -e "Currently processing:";
find . -type f -name \*.SAM -print | xargs -I{} -n 1 -P ${config[CORE]} perl "$PERLDIR" {} ${config[GENOME_NAME]} ${config[SPACE_SAVER]} ${config[READ1_EXT]} ${config[READ2_EXT]} ${config[LIBRARY_TYPE]} ${config[FASTA_INDEX]} ${config[MAPQ_FILTER]} ${config[REPEAT_LIST]}

#my $fileCheck = scalar(@files);
#print "\nFailed to detect any .SAM files within the current directory.\n\n" if $fileCheck == 0;
#exit if $fileCheck == 0;	#Exit script if no .SAM files are found

#############################
# Histogram Maps & Statistics
#############################
echo "-------------------------------------"
echo "Generating Histogram Maps..."
echo -e "-------------------------------------\n"
mkdir -p "$RES/Histograms"
perl "$BASEDIR/Scripts/histogramMap.pl" "-i" "$RES/Data" "-g" ${config[GENOME_NAME]}
PT=$(date +%s)
alignTime=$((AT-ST)); ((sec=alignTime%60, alignTime/=60, min=alignTime%60, hrs=alignTime/60))
processTime=$((PT-AT)); ((secB=processTime%60, processTime/=60, minB=processTime%60, hrsB=processTime/60))
runtime=$(printf "%d:%02d:%02d" $hrs $min $sec)
runtimeB=$(printf "%d:%02d:%02d" $hrsB $minB $secB)
echo "-----------------------------";
echo -e "Completed"
echo -e "-----------------------------";
echo -e "Alignment Runtime: $runtime";
echo -e "Processing Runtime: $runtimeB";
echo -e "-----------------------------\n";

###########################
# Directory Organisation
###########################
if [ ${config[SPACE_SAVER]} = N ]; then
   mkdir -p "$RES/SAM"
   mv *.SAM ./SAM
elif [ ${config[SPACE_SAVER]} = Y ]; then
   rm -rf $RES/Data
fi
