#!/usr/bin/env bash
#Package Version: 3.0

######################################################################################################
# Author(s): T.J.Cooper
# Updated: 23/11/2017
# Automates batch-processing of FASTQ/SAM files for genome-wide mapping of terminal molecules
######################################################################################################

function usage {
   echo -e "\nUsage: termMapper -i [INPUT FOLDER] -c [CONFIGURATION FILE]\n";
   echo -e "-i INPUT: Input data folder containing paired-end FASTQ files\n"
   echo -e "-c CONFIG: Configuration file specifying user-parameters (termMapper.config)\n"; exit 1;
}

# Dependencies
command -v perl >/dev/null 2>&1 || { echo >&2 "Error: Perl installation not found."; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "Error: Bowtie2 installation not found."; exit 1; }
[ "${BASH_VERSINFO:-0}" -lt 4 ] && { echo "Error: Bash (v4.0+) installation not found."; exit 1; }
BOWTIE="$(bowtie2 --version)"
BTVER=$(echo "$BOWTIE" | awk '/version/ {print $NF}')

# Command-line Options
CONF=""
INPUT=""
while getopts ":c:i:" FLAG; do
   case $FLAG in
      i) INPUT=$OPTARG;;
      c) CONF=$OPTARG;;
      \?) echo -e "\nInvalid option: -$OPTARG"
         usage;;
      :) echo -e "\nOption -$OPTARG requires an argument."
         usage;;
   esac
done
if [ "$#" -eq 0 ]; then
    usage
fi

# Read Config File
BASEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd "$INPUT" || exit
mkdir -p "$INPUT/Logs" || exit
exec > >(tee "$INPUT/Logs/System.log")
shopt -s extglob
declare -A config=()
while IFS='=' read -r k v; do
   [[ $v ]] || continue
   [[ $k = "#"* ]] && continue
   k=${k%+([[:space:]])}
   v=${v#+([[:space:]])}
   config[$k]=$v
done < "$CONF"
if [ "${config[LIBRARY_TYPE]}" = "SINGLE" ]; then
   PERLDIR=$BASEDIR/Scripts/ReadProcessing.pl
elif [ "${config[LIBRARY_TYPE]}" = "OLIGO" ]; then
   PERLDIR=$BASEDIR/Scripts/DualEndExtractOligo.pl
elif [ "${config[LIBRARY_TYPE]}" = "DOUBLE" ]; then
   PERLDIR=$BASEDIR/Scripts/ReadProcessing.pl
fi
if [[ ! -f $PERLDIR || ! -d ${config[GENOME_DIR]} || -z ${config[GENOME_NAME]} || -z ${config[READ1_EXT]}
|| -z ${config[READ2_EXT]} || -z ${config[TRIM]} || -z ${config[TRIM_LEN]} || -z ${config[CORE]} || -z ${config[SPACE_SAVER]} ]]; then
    echo -e "\nError (Config File): User parameters or script files are missing and/or incorrectly specified.\n"; exit -1
fi
if [ ${config[TRIM]} = Y ] && [ ${config[LIBRARY_TYPE]} = OLIGO ]; then
   echo -e "\nError (Config File): Trimmed alignment is not available in OLIGO mode.\n"; exit 1
fi

# Input Files
declare -A STR
shopt -s nullglob
for file in *${config[READ1_EXT]}.*; do
   STR["${file%*${config[READ1_EXT]}.*}"]=1
done
if [ ${#STR[@]} = 0 ]; then
   echo -e "\nNo valid FASTQ files were found."
   usage
fi

# Global Alignment
S=$(date +%s)
echo "--------------------------------------------"
echo "FASTQ Alignment (Global --end-to-end)"
echo "--------------------------------------------"
echo "Currently aligning:"
for i in "${!STR[@]}"; do
   printf "%s\n%s\n%s\n%s\n" "-----------------------" "Bowtie2 Version: $BTVER" "Reference Genome: ${config[GENOME_NAME]}" "Sample: $i" > "$INPUT/Logs/Log.$i.txt"
   echo ${i}
   printf "%s\n%s\n%s\n%s\n" "---------------------------------------------------------------------" "GLOBAL" "Parameters: ${config[GLOBAL_OPTIONS]}" "---------------------------------------------------------------------" >> "$INPUT/Logs/Log.$i.txt"
   bowtie2 --end-to-end ${config[GLOBAL_OPTIONS]} -p ${config[CORE]} -x ${config[GENOME_DIR]}/${config[GENOME_NAME]} -1 $i${config[READ1_EXT]}.fastq -2 $i${config[READ2_EXT]}.fastq -S $i"_Global".SAM 2>> "$INPUT/Logs/Log.$i.txt"
   if [ ${config[SPACE_SAVER]} = Y ]; then
      rm $i${config[READ1_EXT]}.fastq
      rm $i${config[READ2_EXT]}.fastq
   fi
done
perl "$PERLDIR" "Global.SAM" ${config[SPACE_SAVER]} ${config[GENOME_NAME]} ${config[TRIM]} ${config[TRIM_LEN]} ${config[READ1_EXT]} ${config[READ2_EXT]} ${config[LIBRARY_TYPE]} ${config[REPEAT_LIST]}
E=$(date +%s)
i=$((E-S)); ((sec=i%60, i/=60, min=i%60, hrs=i/60))
runtime=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "-----------------------";
echo -e "Completed\nRuntime: $runtime";
echo -e "-----------------------\n";

# Trimmed Alignment
if [ ${config[TRIM]} = Y ]; then
   echo "------------------------------------------"
   echo "FASTQ Alignment (Trimmed --local)"
   echo "------------------------------------------"
   echo "Currently realigning:"
   for n in "${!STR[@]}"; do
      echo ${n}
      printf "\n%s\n%s\n%s\n%s\n" "---------------------------------------------------------------------" "TRIMMED" "Parameters: ${config[GLOBAL_OPTIONS]}" "---------------------------------------------------------------------" >> "$INPUT/Logs/Log.${n}.txt"
      bowtie2 --local ${config[LOCAL_OPTIONS]} -p ${config[CORE]} -x ${config[GENOME_DIR]}/${config[GENOME_NAME]} -1 ${n}${config[READ1_EXT]}"_unmapped_trimmed.fastq" -2 ${n}${config[READ2_EXT]}"_unmapped_trimmed.fastq" -S $n"_Trimmed".SAM 2>> "$INPUT/Logs/Log.${n}.txt"
      if [ ${config[SPACE_SAVER]} = Y ]; then
         rm ${n}${config[READ1_EXT]}"_unmapped_trimmed.fastq"
         rm ${n}${config[READ2_EXT]}"_unmapped_trimmed.fastq"
      fi
   done
   perl "$PERLDIR" "Trimmed.SAM" ${config[SPACE_SAVER]} ${config[GENOME_NAME]}
   ET=$(date +%s)
   it=$((ET-E)); ((sec=it%60, it/=60, min=it%60, hrs=it/60))
   runtime=$(printf "%d:%02d:%02d" $hrs $min $sec)
   echo "-----------------------";
   echo -e "Completed\nRuntime: $runtime";
   echo -e "-----------------------\n";
fi

# Directory Organisation
if [ ${config[SPACE_SAVER]} = N ]; then
   mkdir -p "$INPUT/FASTQ"
   mkdir -p "$INPUT/SAM"
   mv *.fastq ./FASTQ
   mv *.SAM ./SAM
fi

# Histogram Maps & Statistics
echo "------------------------------------------"
echo "Generating Histogram Maps..."
echo -e "------------------------------------------\n"
cd "$INPUT" || exit
mkdir -p "$INPUT/Histograms"
perl "$BASEDIR/Scripts/HistogramMap.pl" "-i" "$INPUT/Data" "-g" ${config[GENOME_NAME]}
