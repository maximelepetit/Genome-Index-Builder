#!/bin/bash
set -e

# Function to replace spaces with underscores and capitalize the first letter of each word
transform_and_capitalize() {
    echo "${1// /_}" | awk '{for(i=1;i<=NF;i++){$i=toupper(substr($i,1,1)) tolower(substr($i,2));}print}'
}

# Logging function for standardized message format
log() {
    echo "$(date) - $1"
}

usage() {
    log " Usage: $0 -PathOutputReference <dir> -speciesNames <species1> [<species2>] -condaEnv <name> [-threads <num> -PathGenFastaFile <path> [<path>] -PathGtfFile <path> [<path>]]
    
    Required arguments:
    -PathOutputReference : Directory where reference will be created.
    -speciesNames : Organism type of sample, usually referring to species.
    -condaEnv : Conda environment with Split-Pipe installed.

    Optional arguments:
    -threads : Number of threads to use.
    -PathGenFastaFile : Path to genome FASTA file(s). [1 or 2 paths]
    -PathGtfFile : Path to GTF annotation file(s). [1 or 2 paths]
    "
    exit 1
}
sort_array() {
    local array=("$@")
    local sorted

    if ! sorted=$(printf "%s\n" "${array[@]}" | sort); then
        printf "Error: Sorting failed.\n" >&2
        return 1
    fi

    printf "%s\n" "$sorted"
}  
activate_conda() {
    # Ensure the conda initialization script is sourced
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        . "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        . "/opt/conda/etc/profile.d/conda.sh"
    else
        log " - Cannot find conda initialization script."
        exit 1
    fi
    conda activate $condaEnv || { log " - Failed to activate conda environment: $condaEnv"; exit 1; }
}
deactivate_conda() {
    # Ensure the conda initialization script is sourced
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        . "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        . "/opt/conda/etc/profile.d/conda.sh"
    else
        log " - Cannot find conda initialization script."
        exit 1
    fi

    conda deactivate  || { log " - Failed to deactivate conda environment: $condaEnv"; exit 1; }
}

# Initialize variables
condaEnv=splitpipe-1.5.0
LC_TIME=fr_FR.UTF-8
speciesNames=()
fastaFiles=()
gtfFiles=()



while [[ -n "$1" ]]; do
    case "$1" in
        -condaEnv) condaEnv="$2"; shift 2 ;;
        -speciesNames) 
            shift
            while [[ -n "$1" && ! "$1" =~ ^- ]]; do
                speciesNames+=("$1")
                shift
            done
            ;;
        -PathOutputReference) PathOutputReference="$2"; shift 2 ;;
        -threads) threads="$2"; shift 2 ;;
        -PathGenFastaFile) 
            shift
            while [[ -n "$1" && ! "$1" =~ ^- ]]; do
                fastaFiles+=("$1")
                shift
            done
            ;;
        -PathGtfFile) 
            shift
            while [[ -n "$1" && ! "$1" =~ ^- ]]; do
                gtfFiles+=("$1")
                shift
            done
            ;;
        *) log " - Unknown option: $1"; usage ;;
    esac
done



# Check for mandatory arguments
if [[ -z "$PathOutputReference" || -z "$speciesNames" || -z "$condaEnv" ]]; then
    log " - Error: Missing required arguments."
    usage
fi

# Check the number of genome names and species names
if [[ ${#speciesNames[@]} -gt 2 ]]; then
    log " - Too many species names provided."
    exit 1
fi


# Ensure the number of FASTA and GTF files provided match and are valid
if [[ ${#fastaFiles[@]} -ne 0 || ${#gtfFiles[@]} -ne 0 ]]; then
    # If only FASTA or only GTF files are provided, but not both, show an error
    if [[ ${#fastaFiles[@]} -eq 0 || ${#gtfFiles[@]} -eq 0 ]]; then
        log " - Error: If providing FASTA or GTF files, you must provide both."
        exit 1
    fi

    # Ensure equal numbers of FASTA, GTF, and species names
    if [[ ${#fastaFiles[@]} -ne ${#gtfFiles[@]} || ${#fastaFiles[@]} -ne ${#speciesNames[@]} || ${#gtfFiles[@]} -ne ${#speciesNames[@]} ]]; then
        log " - Error: The number of FASTA and GTF files must match and correspond to the number of species names."
        exit 1
    fi

    # Ensure no more than 2 FASTA and 2 GTF files
    if [[ ${#fastaFiles[@]} -gt 2 || ${#gtfFiles[@]} -gt 2 ]]; then
        log " - Error: You can provide a maximum of 2 FASTA and 2 GTF files."
        exit 1
    fi
fi






# Check and transform species names if provided
for i in "${!speciesNames[@]}"; do
    if [[ ! "${speciesNames[$i]}" =~ ^[a-zA-Z_[:space:]]+$ ]]; then
        log " - Species name must contain only letters or underscores (Mus_musculus, mus musculus, mus_musculus)"
        exit 1
    else
        speciesNames[$i]=$(transform_and_capitalize "${speciesNames[$i]}")
        log " - Species: ${speciesNames[$i]}"
    fi
done

if [[ -z  $threads ]] ; then
    threads=$(( $(nproc) / 2 ))
else 
    # Check if threads is a valid positive integer
    if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
        log " - Threads '$threads' is not a valid positive integer"
        exit 1
    fi
fi 

log " - Threads: $threads"

# Directories
dateSuffix=$(date '+%Y_%m_%d')
referenceDir="${PathOutputReference}/${speciesNames[0]}${speciesNames[1]:+_${speciesNames[1]}}/$dateSuffix"


if [ -d "${referenceDir}" ]; then
    log " - Directory ${referenceDir} already exist. Remove ${referenceDir}"
    exit 1
else
    log " - Directory ${referenceDir} does not exist. Creating directory."
    mkdir -p "$referenceDir"
    log " - Created new subdirectory: $referenceDir"
fi

# Use referenceDir for further operations
log " - Using $referenceDir for your operations."




# Create necessary directories
mkdir -p "$referenceDir/genome/dna"  "$referenceDir/genes" 

cd "$referenceDir" || {
    log " - Failed to change directory to $referenceDir"
    exit 1
}






# Logic for handling species mix with provided files or single species
if [[ ${#fastaFiles[@]} -ne 0 && ${#gtfFiles[@]} -ne 0 ]]; then
    for fasta in "${fastaFiles[@]}"; do
        cp "$fasta" "$referenceDir/genome/dna"
        if [[ "$fasta" != *.gz ]]; then
            gzip -f "$referenceDir/genome/dna/$(basename "$fasta")"
        fi
    done
    for gtf in "${gtfFiles[@]}"; do
        cp "$gtf" "$referenceDir/genes/"
        if [[ "$gtf" != *.gz ]]; then
            gzip -f "$referenceDir/genes/$(basename "$gtf")"
        fi
    done

else
    for species in "${speciesNames[@]}"; do

        ensembl_species="${species,,}"
        request_gtf_fd="ftp://ftp.ensembl.org/pub/current_gtf/${ensembl_species}/"
        request_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/dna/"
        
                
        log " - Trying to download primary assembly file for $species from Ensembl..."
        genome_filename="*.dna.primary_assembly.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna" -A "$genome_filename" "$request_fasta_fd" || {
        log " - Primary assembly file not found for $species. Trying to download toplevel file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna" -A "$genome_filename" "$request_fasta_fd" || {
        log " - Failed to download DNA file for $species from Ensembl."
        exit 1
        }
        }
        log " - Trying to download GTF file for $species from Ensembl..."
        gtf_filename="*[0-9].gtf.gz"
        wget -r -np -nd -q -P "$referenceDir/genes" -A "${gtf_filename}" "$request_gtf_fd" || {
        log " - Failed to download GTF file for $species from Ensembl."
        exit 1
        }
    done

fi



 
genomeFastaFiles=($(find "${referenceDir}/genome/dna/" -type f -name "*.fa.gz"))
annotationGTFFiles=($(find "${referenceDir}/genes/" -type f -name "*.gtf.gz"))



speciesNames=($(sort_array "${speciesNames[@]}"))

annotationGTFFiles=($(sort_array "${annotationGTFFiles[@]}"))

genomeFastaFiles=($(sort_array "${genomeFastaFiles[@]}"))



log " - Creating Split-pipe index..."
activate_conda


/usr/bin/time split-pipe \
    --mode mkref \
    --nthreads "$threads" \
    --genome_name "${speciesNames[@]}" \
    --fasta "${genomeFastaFiles[@]}" \
    --genes "${annotationGTFFiles[@]}" \
    --output_dir "${referenceDir}/makeRef"

deactivate_conda

log " - Split-pipe Index complete."
rm -rf "${referenceDir}/genome/" "$referenceDir/genes/"