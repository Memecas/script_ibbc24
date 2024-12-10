#!/bin/bash

# Color definitions
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
RED="\033[0;31m"
BLUE="\033[0;34m"
RESET="\033[0m" # Reset color

# Default parameters for Fastp
DEFAULT_QUALITY=20
DEFAULT_MIN_LENGTH=50
DEFAULT_ADAPTER_DETECTION="enabled"
DEFAULT_COMPLEXITY_FILTER="disabled"
DEFAULT_THREADS=8

# Display parameters and prompt for confirmation
echo -e "${BLUE}The script will run with the following default parameters for fastp:${RESET}"
echo -e "${YELLOW}1. Quality threshold (-q): $DEFAULT_QUALITY (Bases with quality below this value will be trimmed)"
echo -e "2. Minimum read length (-l): $DEFAULT_MIN_LENGTH (Reads shorter than this length will be discarded)"
echo -e "3. Adapter detection (--detect_adapter_for_pe): $DEFAULT_ADAPTER_DETECTION (Automatically detect adapters for paired-end data)"
echo -e "4. Complexity filter (--low_complexity_filter): $DEFAULT_COMPLEXITY_FILTER (Remove low-complexity reads)"
echo -e "5. Number of threads (--thread): $DEFAULT_THREADS (Number of threads for parallel processing)${RESET}"
echo ""
echo -e "${BLUE}Do you want to proceed with these parameters? (yes/no)${RESET}"

read -p "Enter your choice: " confirm

if [[ "$confirm" != "yes" ]]; then
    echo -e "${YELLOW}You have chosen to modify the parameters. Please provide new values for each parameter as prompted below.${RESET}"

    # Prompt for each parameter
    read -p "1. Quality threshold (-q) [Current: $DEFAULT_QUALITY]: " quality
    [[ -z "$quality" ]] && quality=$DEFAULT_QUALITY

    read -p "2. Minimum read length (-l) [Current: $DEFAULT_MIN_LENGTH]: " min_length
    [[ -z "$min_length" ]] && min_length=$DEFAULT_MIN_LENGTH

    read -p "3. Enable adapter detection (--detect_adapter_for_pe)? (yes/no) [Current: $DEFAULT_ADAPTER_DETECTION]: " adapter_detection
    if [[ "$adapter_detection" == "no" ]]; then
        adapter_detection="disabled"
    else
        adapter_detection="enabled"
    fi

    read -p "4. Enable complexity filter (--low_complexity_filter)? (yes/no) [Current: $DEFAULT_COMPLEXITY_FILTER]: " complexity_filter
    if [[ "$complexity_filter" == "yes" ]]; then
        complexity_filter="enabled"
    else
        complexity_filter="disabled"
    fi

    read -p "5. Number of threads (--thread) [Current: $DEFAULT_THREADS]: " threads
    [[ -z "$threads" ]] && threads=$DEFAULT_THREADS

    # Update parameters with user inputs
    DEFAULT_QUALITY=$quality
    DEFAULT_MIN_LENGTH=$min_length
    DEFAULT_ADAPTER_DETECTION=$adapter_detection
    DEFAULT_COMPLEXITY_FILTER=$complexity_filter
    DEFAULT_THREADS=$threads

    echo ""
    echo -e "${YELLOW}Updated parameters:${RESET}"
    echo -e "1. Quality threshold (-q): $DEFAULT_QUALITY"
    echo -e "2. Minimum read length (-l): $DEFAULT_MIN_LENGTH"
    echo -e "3. Adapter detection (--detect_adapter_for_pe): $DEFAULT_ADAPTER_DETECTION"
    echo -e "4. Complexity filter (--low_complexity_filter): $DEFAULT_COMPLEXITY_FILTER"
    echo -e "5. Number of threads (--thread): $DEFAULT_THREADS"
fi

cat << "EOF"
   _____ ____  _   _ _____          
  / ____/ __ \| \ | |  __ \   /\    
 | |   | |  | |  \| | |  | | /  \   
 | |   | |  | | . ` | |  | |/ /\ \  
 | |___| |__| | |\  | |__| / ____ \ 
  \_____\____/|_| \_|_____/_/    \_\
                                    
                                         
                              
EOF

# Check if the correct conda environment is activated
EXPECTED_ENV="tools_qc"
if [[ "$CONDA_DEFAULT_ENV" != "$EXPECTED_ENV" ]]; then
    echo -e "${RED}Warning: The expected conda environment '$EXPECTED_ENV' is not active.${RESET}"

    # Try to activate the environment
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$EXPECTED_ENV"

    # Check if activation was successful
    if [[ "$CONDA_DEFAULT_ENV" == "$EXPECTED_ENV" ]]; then
        echo -e "${GREEN}Successfully activated the conda environment: $CONDA_DEFAULT_ENV${RESET}"
    else
        echo -e "${RED}Error: Failed to activate the conda environment '$EXPECTED_ENV'. Please activate it manually and re-run the script.${RESET}"
        exit 1
    fi
else
    echo -e "${GREEN}Conda environment detected: $CONDA_DEFAULT_ENV${RESET}"
fi

# Check for dependencies
MISSING_DEPENDENCIES=()
for tool in fastqc fastp multiqc; do
    if ! command -v "$tool" &> /dev/null; then
        MISSING_DEPENDENCIES+=("$tool")
    fi
done

if [[ ${#MISSING_DEPENDENCIES[@]} -ne 0 ]]; then
    echo -e "${RED}Error: The following dependencies are missing in the current environment ($CONDA_DEFAULT_ENV): ${MISSING_DEPENDENCIES[*]}${RESET}"
    echo "Please install them and re-run the script."
    exit 1
fi

# Base directory setup
BASE_DIR=$(pwd)
INPUT_DIR="$BASE_DIR"
RESULTS_DIR="$BASE_DIR/results"
FASTQC_BEFORE_DIR="$RESULTS_DIR/fastqc/before_cleaning"
FASTQC_AFTER_DIR="$RESULTS_DIR/fastqc/after_cleaning"
FASTP_CLEANED_DIR="$RESULTS_DIR/fastp_cleaned" 
MULTIQC_REPORT_DIR="$RESULTS_DIR/multiqc"
REPORTS_DIR="$RESULTS_DIR/reports"

# Function to create directories with error handling
create_directories() {
    echo -e "${BLUE}Creating directories...${RESET}"
    dirs=("$FASTQC_BEFORE_DIR" "$FASTQC_AFTER_DIR" "$FASTP_CLEANED_DIR" "$MULTIQC_REPORT_DIR" "$REPORTS_DIR")

    for dir in "${dirs[@]}"; do
        if mkdir -p "$dir"; then
            echo -e "${GREEN}Created directory: $dir${RESET}"
        else
            echo -e "${RED}Error: Failed to create directory: $dir${RESET}"
            exit 1
        fi
    done
}

#Function to process samples
process_samples() {
    echo -e "${BLUE}Processing samples...${RESET}"
    echo "Sample_Name,Initial_Reads,Processed_Reads,Discarded_Reads" > "$REPORTS_DIR/read_summary.csv"

    # Iterate over input files to identify paired-end or single-end samples
    for file in "$INPUT_DIR"/*_1_*.fastq.gz "$INPUT_DIR"/*R1*.fastq.gz "$INPUT_DIR"/*_1.fastq.gz "$INPUT_DIR"/*R1.fastq; do
        if [[ -f "$file" ]]; then
            sample_name=$(basename "$file" | cut -d'_' -f1,2)
            r1="$file"
            r2="${r1/_1_/_2_}"  # Replace _1_ with _2_
            r2="${r2/_1.fastq/_2.fastq}"  # Adjust for .fastq naming convention
            r2="${r2/R1/R2}"  # Handle R1 -> R2

            # Check if R2 exists or handle as single-end
            if [[ ! -f "$r2" ]]; then
                echo -e "${YELLOW}Warning: Paired-end file not found for $r1. Processing as single-end.${RESET}"
                r2=""
            fi

cat << "EOF"
   __          _    ____   _____ 
  / _|        | |  / __ \ / ____|
 | |_ __ _ ___| |_| |  | | |     
 |  _/ _` / __| __| |  | | |     
 | || (_| \__ \ |_| |__| | |____ 
 |_| \__,_|___/\__|\___\_\\_____|
                                 
                                                                                   
                              
EOF

            # Check if FastQC has already been run
            fastqc_output_dir="$FASTQC_BEFORE_DIR/$sample_name"
            if [[ -d "$fastqc_output_dir" && "$(ls -A "$fastqc_output_dir")" ]]; then
                echo -e "${YELLOW}FastQC reports already exist for $sample_name. Skipping FastQC before cleaning.${RESET}"
            else
                echo -e "${BLUE}Running FastQC for raw files: $sample_name...${RESET}"
                mkdir -p "$fastqc_output_dir"
                fastqc -o "$fastqc_output_dir" "$r1" "$r2"
            fi

cat << "EOF"
   __          _   _____  
  / _|        | | |  __ \ 
 | |_ __ _ ___| |_| |__) |
 |  _/ _` / __| __|  ___/ 
 | || (_| \__ \ |_| |     
 |_| \__,_|___/\__|_|     
                          
                           
                              
EOF

            # Run Fastp
            echo -e "${BLUE}Running Fastp for $sample_name...${RESET}"
            fastp_output_prefix="$FASTP_CLEANED_DIR/$sample_name"
            fastp_log="$REPORTS_DIR/${sample_name}_fastp.log"

            fastp --in1 "$r1" --in2 "$r2" --out1 "${fastp_output_prefix}_cleaned_1.fastq.gz" \
                --out2 "${fastp_output_prefix}_cleaned_2.fastq.gz" \
                --thread $DEFAULT_THREADS --qualified_quality_phred $DEFAULT_QUALITY \
                --length_required $DEFAULT_MIN_LENGTH \
                $( [[ "$DEFAULT_ADAPTER_DETECTION" == "enabled" ]] && echo "--detect_adapter_for_pe" ) \
                $( [[ "$DEFAULT_COMPLEXITY_FILTER" == "enabled" ]] && echo "--low_complexity_filter" ) \
                --html "$fastp_output_prefix.html" --json "$fastp_output_prefix.json" \
                &> "$fastp_log"

            # Add Read Summary
            initial_reads=$(grep "total reads" "$fastp_log" | head -1 | awk '{print $NF}')
            processed_reads=$(grep "total reads" "$fastp_log" | tail -1 | awk '{print $NF}')
            discarded_reads=$((initial_reads - processed_reads))
            echo "$sample_name,$initial_reads,$processed_reads,$discarded_reads" >> "$REPORTS_DIR/read_summary.csv"

cat << "EOF"
   __          _    ____   _____ 
  / _|        | |  / __ \ / ____|
 | |_ __ _ ___| |_| |  | | |     
 |  _/ _` / __| __| |  | | |     
 | || (_| \__ \ |_| |__| | |____ 
 |_| \__,_|___/\__|\___\_\\_____|
                                 
                                                                          
                              
EOF

            # FastQC after cleaning
            fastqc_after_output_dir="$FASTQC_AFTER_DIR/$sample_name"
            if [[ -d "$fastqc_after_output_dir" && "$(ls -A "$fastqc_after_output_dir")" ]]; then
                echo -e "${YELLOW}FastQC reports already exist for cleaned files of $sample_name. Skipping FastQC after cleaning.${RESET}"
            else
                echo -e "${BLUE}Running FastQC for cleaned files: $sample_name...${RESET}"
                mkdir -p "$fastqc_after_output_dir"
                fastqc -o "$fastqc_after_output_dir" "${fastp_output_prefix}_cleaned_1.fastq.gz" "${fastp_output_prefix}_cleaned_2.fastq.gz"
            fi
        fi
    done
}

# Function to generate MultiQC report
generate_multiqc_report() {
    echo -e "${BLUE}Generating MultiQC report...${RESET}"
    multiqc "$RESULTS_DIR/fastqc" -o "$MULTIQC_REPORT_DIR"
}

# Main script execution
create_directories
process_samples
generate_multiqc_report

echo -e "${GREEN}Pipeline completed successfully. Reports and cleaned files are available in the results directory.${RESET}"
