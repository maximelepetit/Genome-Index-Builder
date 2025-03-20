#!/bin/bash
set -e
#bash ./dev/preBuild_Genomev1.1.sh -outDir "/mnt/Data1/genomes/stomics" -species "Gallus_gallus" -sif "/home/labex-cortex/singularity_img/SAW_7.1.sif" -threads 20 -fasta "/mnt/Data1/genomes/gallus_gallus/GRCg6a/gallus_GRCg6a/GalGal_GRCg6a.fa.gz" -gtf "/mnt/Data1/genomes/gallus_gallus/GRCg6a/gallus_GRCg6a/GalGal_GRCg6a.gtf.gz" 
# Function to replace spaces with underscores and capitalize the first letter of each word
transform_and_capitalize() {
    echo "${1// /_}" | awk '{for(i=1;i<=NF;i++){$i=toupper(substr($i,1,1)) tolower(substr($i,2));}print}'
}

usage() {
    echo "$(date) usage: $0 -PathOutputReference <path> -species <name> [-threads <num> -PathGenFastaFile <path> -PathGtfFile <path>]
    -PathOutputReference : Directory where reference will be created (required)
    -species : Organism type of sample, usually referring to species. (required)
    -threads : the number of threads to be used in running the pipeline (default: 4) (optional)
    -PathGenFastaFile : path to custom .fa file (optional)
    -PathGtfFile : path to custom .gtf file (optional)
    "
    exit 1
}
# Initialize variables
threads=4
LC_TIME=fr_FR.UTF-8



while [[ -n "$1" ]]; do
    case "$1" in
        -PathOutputReference) PathOutputReference="$2"; shift ;;
        -species) species="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathGtfFile) PathGtfFile="$2"; shift ;;
        *) echo "$(date) Unknown option: $1"; usage ;;
    esac
    shift
done

# Check for mandatory arguments
if [[ -z "$PathOutputReference" || -z "$species" ]]; then
    usage
fi


# Check speciesName or custom files
if [[ -z "$species" && (-z "$PathGenFastaFile" || -z "$PathGtfFile") ]]; then
    echo "$(date) - Either species name or both fasta and gtf files must be provided."
    exit 1
fi

# If fasta and gtf are provided, validate them
if [[ -n "$PathGenFastaFile" && -n "$PathGtfFile" ]]; then
    if [[ ! -f "$PathGenFastaFile" ]]; then
        echo "$(date) - Provided DNA fasta file does not exist: $PathGenFastaFile"
        exit 1
    fi
    if [[ ! -f "$PathGtfFile" ]]; then
        echo "$(date) - Provided GTF file does not exist: $PathGtfFile"
        exit 1
    fi
else
    # If only one of the custom files is provided, exit with error
    if [[ -n "$PathGenFastaFile" || -n "$PathGtfFile" ]]; then
        echo "$(date) - Both DNA fasta and GTF files must be provided together."
        exit 1
    fi
fi

# Check and transform species name if provided
if [[ -n "$species" ]]; then
    if [[ ! "$species" =~ ^[a-zA-Z_[:space:]]+$ ]]; then
        echo "$(date) - Species name must contain only letters or underscores (Mus musculus, mus musculus, mus_musculus)"
        exit 1
    else
        species=$(transform_and_capitalize "$species")
        echo "$(date) - Species: $species"
    fi
fi


# Check if threads is a valid positive integer
if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
    echo "$(date) - Threads '$threads' is not a valid positive integer"
    exit 1
fi

echo "$(date) - Threads: $threads"

# Directories
dateSuffix=$(date '+%Y_%m_%d')
referenceDir="$PathOutputReference/$species/$dateSuffix"

if [ -d "$referenceDir" ]; then
    echo "Directory $referenceDir already exist. Remove $referenceDir please"
    exit 1
else
    echo "Directory $referenceDir does not exist. Creating directory."
    mkdir -p "$referenceDir"
    echo "Created new subdirectory: $referenceDir"
fi

# Use referenceDir for further operations
echo "Using $referenceDir for your operations."


# Create necessary directories
mkdir -p "$referenceDir/genome" "$referenceDir/genes"

cd "$referenceDir" || {
    echo "$(date) - Failed to change directory to $referenceDir"
    exit 1
}

# Logic for custom files
if [[ -n "$PathGenFastaFile" && -n "$PathGtfFile" ]]; then
    cp "$PathGenFastaFile" "$referenceDir/genome/"
    cp "$PathGtfFile" "$referenceDir/genes/"
    
    if [[ "$PathGenFastaFile" == *.gz ]]; then
        gunzip -f "$referenceDir/genome/$(basename "$PathGenFastaFile")"
        genomeFastaFiles=$(find "${referenceDir}/genome/" -type f -name "$(basename "$PathGenFastaFile" .gz)" | head -n 1)
    else
        genomeFastaFiles=$(find "${referenceDir}/genome/" -type f -name "$(basename "$PathGenFastaFile")" | head -n 1)
    fi

    if [[ "$PathGtfFile" == *.gz ]]; then
        gunzip -f "$referenceDir/genes/$(basename "$PathGtfFile")"
        pathFileGTF=$(find "${referenceDir}/genes/" -type f -name "$(basename "$PathGtfFile" .gz)" | head -n 1)
    else
        pathFileGTF=$(find "${referenceDir}/genes/" -type f -name "$(basename "$PathGtfFile")" | head -n 1)
    fi

    pathFileCheckGTF="${referenceDir}/genes/$(basename "$pathFileGTF" .gtf).check.gtf"

    echo "$(date) - Checking GTF file format..."

    /usr/bin/time saw checkGTF \
    --input-gtf="${pathFileGTF}" \
    --output-gtf="${pathFileCheckGTF}"

    echo "$(date) - GTF file check complete."

    

    echo "$(date) - Creating genome..."
    /usr/bin/time saw makeRef \
        --mode=STAR \
        --genome="${referenceDir}/makeRef" \
        --fasta="${genomeFastaFiles}" \
        --gtf="${pathFileGTF}" \
        --threads-num="$threads"


else
    ensembl_species="${species,,}"

    request_gtf_fd="ftp://ftp.ensembl.org/pub/current_gtf/${ensembl_species}/"
    request_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/dna/"

    gtf_filename="*[0-9].gtf.gz"
    genome_filename="*.dna.primary_assembly.fa.gz"
                
    echo "$(date) - Trying to download primary assembly file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genome" \
        -A "$genome_filename" \
        "$request_fasta_fd" || {
        echo "$(date) - Primary assembly file not found. Trying to download toplevel file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome" \
            -A "$genome_filename" \
            "$request_fasta_fd" || {
            echo "$(date) - Failed to download DNA file from Ensembl."
            exit 1
        }
    }

    echo "$(date) - Trying to download GTF file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genes" -A "${gtf_filename}" "$request_gtf_fd" || {
        echo "$(date) - Failed to download GTF file from Ensembl."
        exit 1
    }

    echo "$(date) - Unzip downloaded files"
    gunzip -f "$referenceDir/genome/"*.fa.gz || {
        echo "$(date) - Failed to unzip DNA files."
        exit 1
    }

    gunzip -f "$referenceDir/genes/"*.gtf.gz || {
        echo "$(date) - Failed to unzip GTF files."
        exit 1
    }

    echo "$(date) - The FASTA and GTF files have been downloaded successfully"

    pathFileGTF=$(find "${referenceDir}/genes/" -type f -regex "${referenceDir}/genes/"*.gtf | head -n 1)
   
    pathFileCheckGTF="${referenceDir}/genes/$(basename "$pathFileGTF" .gtf).check.gtf"

    genomeFastaFiles=$(find "${referenceDir}/genome/" -type f -regex "${referenceDir}/genome/${species}.${assemblyVersion}.*.dna.\(primary_assembly\|toplevel\).fa" | head -n 1)


    echo "$(date) - Checking GTF file format..."
    /usr/bin/time saw checkGTF \
    --input-gtf="${pathFileGTF}" \
    --output-gtf="${pathFileCheckGTF}"
    echo "$(date) - GTF file check complete."


    echo "$(date) - Creating genome..."
    /usr/bin/time saw makeRef \
        --mode=STAR \
        --genome="${referenceDir}/makeRef" \
        --fasta="${genomeFastaFiles}" \
        --gtf="${pathFileCheckGTF}" \
        --threads-num="$threads"
    
fi
