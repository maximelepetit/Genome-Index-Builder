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
    log " usage: $0 -PathOutputReference <path> -speciesNames <name> [-threads <num> -PathGenFastaFile <path> -PathGtfFile <path> -checkGTF  ]
    
Required arguments:
    -PathOutputReference   : Directory where reference will be created.
    -speciesNames          : Organism type of sample, usually referring to species. 

Optional arguments:
    -threads               : Number of threads to use.
    -memory                : Maximum memory [Gb] allocated for the pipeline . (default: free / 2)
    -PathGenFastaFile      : Path to genome FASTA file.
    -PathGtfFile           : Path to GTF annotation file.
    -checkGTF              : Flag to check GTF annotation file format.
    "
    exit 1
}
# Initialize variables
threads=4
LC_TIME=fr_FR.UTF-8



while [[ -n "$1" ]]; do
    case "$1" in
        -PathOutputReference) PathOutputReference="$2"; shift ;;
        -speciesNames) speciesNames="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -memory) memory="$2"; shift ;;
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathGtfFile) PathGtfFile="$2"; shift ;;
        -checkGTF) checkGTF=1 ;;  
        *) log " - Unknown option: $1"; usage ;;
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
        log " - Provided DNA fasta file does not exist: ${PathGenFastaFile}"
        exit 1
    fi
    if [[ ! -f "$PathGtfFile" ]]; then
        log " - Provided GTF file does not exist: ${PathGtfFile}"
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
        log " - Species: ${speciesNames}"
    fi
fi


if [[ -z  $threads ]] ; then

    threads=$(( $(nproc) / 2 ))

else 
    # Check if threads is a valid positive integer
    if [[ ! "$threads" =~ ^[0-9]+$ ]]; then
        log " - Threads '${threads}' is not a valid positive integer"
        exit 1
    fi
fi 

log " - Threads: $threads"

if [[ -z  $memory ]] ; then
    memory=$(( $(free | grep Mem | awk '{print $7}') / 1000000 / 4 ))
else 
    # Check if threads is a valid positive integer
   if ! [[ "$memory" =~ ^[1-9][0-9]*$ ]]; then
    log "Error: Memory must be a positive integer."
    exit 1
    fi
fi

log " - Memory used: ${memory} Gb"


# Directories
dateSuffix=$(date '+%Y_%m_%d')
referenceDir="$PathOutputReference/$speciesNames/$dateSuffix"

if [ -d "$referenceDir" ]; then
    log " - Directory ${referenceDir} already exist. Remove ${referenceDir} please"
    exit 1
else
    log " - Directory ${referenceDir} does not exist. Creating directory."
    mkdir -p "$referenceDir"
    log " - Created new subdirectory: ${referenceDir}"
fi

# Use referenceDir for further operations
log " - Using ${referenceDir} for your operations."


# Create necessary directories
mkdir -p "${referenceDir}/genome/dna" "${referenceDir}/genes"



# Logic for custom files
if [[ -n "$PathGenFastaFile" && -n "$PathGtfFile" ]]; then
    cp "$PathGenFastaFile" "$referenceDir/genome/dna"
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
    wget -r -np -nd -q -P "$referenceDir/genome/dna" \
        -A "$genome_filename" \
        "$request_fasta_fd" || {
        log " - Primary assembly file not found. Trying to download toplevel file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna" \
            -A "$genome_filename" \
            "$request_fasta_fd" || {
            log " - Failed to download DNA file from Ensembl."
            exit 1
        }
    }

    log " - Trying to download GTF file from Ensembl..."
    wget -r -np -nd -q -P "${referenceDir}/genes" -A "${gtf_filename}" "$request_gtf_fd" || {
        log " - Failed to download GTF file from Ensembl."
        exit 1
    }

    log " - Unzip downloaded files"
    gunzip -f "${referenceDir}/genome/dna/"*.fa.gz || {
        log " - Failed to unzip DNA files."
        exit 1
    }

    gunzip -f "${referenceDir}/genes/"*.gtf.gz || {
        log " - Failed to unzip GTF files."
        exit 1
    }

    log " - The FASTA and GTF files have been downloaded successfully"
    
fi

annotationGTFFiles=($(find "${referenceDir}/genes/" -type f -name "*.gtf"))
genomeFastaFiles=($(find "${referenceDir}/genome/dna/" -type f -name "*.fa"))

if [[ "$checkGTF" -eq 1 ]]; then  # Run ONLY if -checkGTF flag is provided

    annotationFilteredGTFFiles="${referenceDir}/genes/$(basename "$annotationGTFFiles" .gtf).check.gtf"

    log " - Running  mkgtf..."


    cellranger mkgtf ${annotationGTFFiles} ${annotationFilteredGTFFiles} \
                    --attribute=gene_biotype:protein_coding \
                    --attribute=gene_biotype:lncRNA \
                    --attribute=gene_biotype:antisense \
                    --attribute=gene_biotype:IG_LV_gene \
                    --attribute=gene_biotype:IG_V_gene \
                    --attribute=gene_biotype:IG_V_pseudogene \
                    --attribute=gene_biotype:IG_D_gene \
                    --attribute=gene_biotype:IG_J_gene \
                    --attribute=gene_biotype:IG_J_pseudogene \
                    --attribute=gene_biotype:IG_C_gene \
                    --attribute=gene_biotype:IG_C_pseudogene \
                    --attribute=gene_biotype:TR_V_gene \
                    --attribute=gene_biotype:TR_V_pseudogene \
                    --attribute=gene_biotype:TR_D_gene \
                    --attribute=gene_biotype:TR_J_gene \
                    --attribute=gene_biotype:TR_J_pseudogene \
                    --attribute=gene_biotype:TR_C_gene


    log " - mkgtf complete."
    annotationGTFFiles=$annotationFilteredGTFFiles
fi



log " - Creating Cellranger index..."
log " - Running  mkref..."
cellranger mkref \
    --output-dir "${referenceDir}/makeRef" \
    --genome ${speciesNames} \
    --fasta ${genomeFastaFiles} \
    --genes ${annotationGTFFiles} \
    --nthreads ${threads} \
    --memgb ${memory}

log " - Cellranger Index complete."
rm -rf "${referenceDir}/genome/" "${referenceDir}/genes/"
