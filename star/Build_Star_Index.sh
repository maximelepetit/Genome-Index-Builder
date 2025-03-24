#!/bin/bash
set -e
#bash ./dev/preBuild_Genomev1.1.sh -outDir "/mnt/Data1/genomes/stomics" -species "Gallus_gallus" -sif "/home/labex-cortex/singularity_img/SAW_7.1.sif" -threads 20 -fasta "/mnt/Data1/genomes/gallus_gallus/GRCg6a/gallus_GRCg6a/GalGal_GRCg6a.fa.gz" -gtf "/mnt/Data1/genomes/gallus_gallus/GRCg6a/gallus_GRCg6a/GalGal_GRCg6a.gtf.gz" 
# Function to replace spaces with underscores and capitalize the first letter of each word
transform_and_capitalize() {
    echo "${1// /_}" | awk '{for(i=1;i<=NF;i++){$i=toupper(substr($i,1,1)) tolower(substr($i,2));}print}'
}

# Logging function for standardized message format
log() {
    echo "$(date) - $1"
}

usage() {
    log " usage: $0 -PathOutputReference <path> -speciesNames <name> [-threads <num> -PathGenFastaFile <path> -PathGtfFile <path>]

    Required arguments:
    -PathOutputReference : Directory where reference will be created.
    -speciesNames : Organism type of sample, usually referring to species.

    Optional arguments:
    -threads : Number of threads to use.
    -PathGenFastaFile : Path to genome FASTA file.
    -PathGtfFile : Path to GTF annotation file to create tx2gene file.
    "
    exit 1
}
# Initialize variables
LC_TIME=fr_FR.UTF-8

while [[ -n "$1" ]]; do
    case "$1" in
        -PathOutputReference) PathOutputReference="$2"; shift ;;
        -speciesNames) speciesNames="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathGtfFile) PathGtfFile="$2"; shift ;;
        *) log " Unknown option: $1"; usage ;;
    esac
    shift
done

# Check for mandatory arguments
if [[ -z "$PathOutputReference" || -z "$speciesNames" ]]; then
    log " - Error: Missing required arguments."
    usage
fi


# Check speciesName or custom files
if [[ -z "$speciesNames" && (-z "$PathGenFastaFile" || -z "$PathGtfFile") ]]; then
    log " - Either species name or both fasta and gtf files must be provided."
    exit 1
fi

# If fasta and gtf are provided, validate them
if [[ -n "$PathGenFastaFile" && -n "$PathGtfFile" ]]; then
    if [[ ! -f "$PathGenFastaFile" ]]; then
        log " - Provided DNA fasta file does not exist: $PathGenFastaFile"
        exit 1
    fi
    if [[ ! -f "$PathGtfFile" ]]; then
        log " - Provided GTF file does not exist: $PathGtfFile"
        exit 1
    fi
else
    # If only one of the custom files is provided, exit with error
    if [[ -n "$PathGenFastaFile" || -n "$PathGtfFile" ]]; then
        log " - Both DNA fasta and GTF files must be provided together."
        exit 1
    fi
fi

# Check and transform species name if provided
if [[ -n "$speciesNames" ]]; then
    if [[ ! "$speciesNames" =~ ^[a-zA-Z_[:space:]]+$ ]]; then
        log " - Species name must contain only letters or underscores (Mus musculus, mus musculus, mus_musculus)"
        exit 1
    else
        speciesNames=$(transform_and_capitalize "$speciesNames")
        log " - Species: $speciesNames"
    fi
fi

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

available_memory=$(( $(free | grep Mem | awk '{print $7}') * 1000 / 4 ))

log " - Memory used: $threads"

# Directories
dateSuffix=$(date '+%Y_%m_%d')
referenceDir="$PathOutputReference/$speciesNames/$dateSuffix"

if [ -d "$referenceDir" ]; then
    log " - Directory $referenceDir already exist. Remove $referenceDir please"
    exit 1
else
    log " - Directory $referenceDir does not exist. Creating directory."
    mkdir -p "$referenceDir"
    log " - Created new subdirectory: $referenceDir"
fi

# Use referenceDir for further operations
log " - Using $referenceDir for your operations."


# Create necessary directories
mkdir -p "$referenceDir/genome/dna"  "$referenceDir/genes"



if [[ -n "$PathGenFastaFile" && -n "$PathGtfFile" ]]; then
    cp "$PathGenFastaFile" "$referenceDir/genome/dna/"
    cp "$PathGtfFile" "$referenceDir/genes/"
    
    if [[ "$PathGenFastaFile" == *.gz ]]; then
        gunzip -f "$referenceDir/genome/dna/$(basename "$PathGenFastaFile")"
    fi

    if [[ "$PathGtfFile" == *.gz ]]; then
        gunzip -f "$referenceDir/genes/$(basename "$PathGtfFile")"
    fi

else
    ensembl_species="${speciesNames,,}"

    request_gtf_fd="ftp://ftp.ensembl.org/pub/current_gtf/${ensembl_species}/"
    request_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/dna/"

    gtf_filename="*[0-9].gtf.gz"
    genome_filename="*.dna.primary_assembly.fa.gz"
                
    log " - Trying to download primary assembly file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genome/dna/" \
        -A "$genome_filename" \
        "$request_fasta_fd" || {
        log " - Primary assembly file not found. Trying to download toplevel file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna/" \
            -A "$genome_filename" \
            "$request_fasta_fd" || {
            log " - Failed to download DNA file from Ensembl."
            exit 1
        }
    }

    log " - Trying to download GTF file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genes" -A "${gtf_filename}" "$request_gtf_fd" || {
        log " - Failed to download GTF file from Ensembl."
        exit 1
    }

    log " - Unzip downloaded files"
    gunzip -f "$referenceDir/genome/dna/"*.fa.gz || {
        log " - Failed to unzip DNA files."
        exit 1
    }

    gunzip -f "$referenceDir/genes/"*.gtf.gz || {
        log " - Failed to unzip GTF files."
        exit 1
    }
   
fi
genomeFastaFiles=($(find "${referenceDir}/genome/dna/" -type f -name "*.fa"))
annotationGTFFiles=($(find "${referenceDir}/genes/" -type f -name "*.gtf"))

log " - Creating Star index..."
/usr/bin/time STAR --runMode genomeGenerate \
    --runThreadN $threads \
    --genomeDir $referenceDir \
    --genomeFastaFiles $genomeFastaFiles \
    --sjdbGTFfile $annotationGTFFiles  \
    --limitGenomeGenerateRAM=$available_memory

log "- Star index complete."
rm -rf "${referenceDir}/genome/" "$referenceDir/genes/"