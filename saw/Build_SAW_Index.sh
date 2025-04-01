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
    log " usage: $0 -PathOutputReference <path> -speciesNames <name> [-threads <num> -PathGenFastaFile <path> -PathGtfFile <path> -rRNA -PathrRNAFastaFile <path> ]
    
    Required arguments:
    -PathOutputReference : Directory where reference will be created.
    -speciesNames : Organism type of sample, usually referring to species. 

    Optional arguments:
    -threads : Number of threads to use.
    -PathGenFastaFile : Path to genome FASTA file.
    -PathGtfFile : Path to GTF annotation file.
    -rRNA : Flag to remove rRNA fragments.
    -PathrRNAFastaFile : Path to GTF annotation file to create tx2gene file.
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
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathGtfFile) PathGtfFile="$2"; shift ;;
        -rRNA) rRNA=1 ;;  
        -PathrRNAFastaFile) PathrRNAFastaFile="$2"; shift ;;
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
mkdir -p "$referenceDir/genome/dna" "$referenceDir/genes"



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

    log " - The FASTA and GTF files have been downloaded successfully"
    
fi


genomeFastaFiles=($(find "${referenceDir}/genome/dna/" -type f -name "*.fa"))
annotationGTFFiles=($(find "${referenceDir}/genes/" -type f -name "*.gtf"))

annotationCheckGTFFiles="${referenceDir}/genes/$(basename "$annotationGTFFiles" .gtf).check.gtf"

log " - Checking GTF file format..."

/usr/bin/time saw checkGTF \
 --input-gtf="${annotationGTFFiles}" \
 --output-gtf="${annotationCheckGTFFiles}"

log " - GTF file check complete."




if [[ "$rRNA" -eq 1 ]]; then  # Run ONLY if -rRNA flag is provided

    mkdir -p "${referenceDir}/genome/rrna/"
    if [[ -f "$PathrRNAFastaFile" && ( "$PathrRNAFastaFile" == *.fa || "$PathrRNAFastaFile" == *.fa.gz ) ]]; then

        log " - Creating tx2gene file from $PathGtfFile"
        cp "$PathrRNAFastaFile" "${referenceDir}/genome/rrna/"
        # Compress if the file is a plain .gtf
        if [[ "$PathrRNAFastaFile" == *.gz ]]; then
            gunzip "${referenceDir}/genome/rrna/$(basename "$PathrRNAFastaFile")"
        fi


    else 
        log " - rRNA FASTA file is invalid or not provided."
        log " - Creating it ..."

        ensembl_species="${speciesNames,,}"  # Ensure lowercase species name
        request_nc_rna_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/ncrna/"
        ncrna_filename="*.ncrna.fa.gz"
        
        wget -r -np -nd -q -P "${referenceDir}/genome/rrna/" -A "$ncrna_filename" "$request_nc_rna_fd" || {
            log " - Failed to download NcRNA FASTA file from Ensembl."
            exit 1

        }
        log " - The NcRNA FASTA files have been downloaded successfully"

        

        PathncRNAFastaFile=$(find "${referenceDir}/genome/rrna/" -type f -name "*.fa.gz" | head -n 1)
        log " - Unzip downloaded files"

        gunzip "${referenceDir}/genome/rrna/$(basename "$PathncRNAFastaFile")"

        PathncRNAFastaFile=$(find "${referenceDir}/genome/rrna/" -type f -name "*.fa" | head -n 1)

        if [[ ! -f "$PathncRNAFastaFile" ]]; then
        echo "ERROR: Input FASTA file $PathncRNAFastaFile not found!"
        exit 1
        fi
        gene_biotype=rRNA

        # Define the output file
        PathrRNAFastaFile="${referenceDir}/genome/rrna/$(basename "$PathncRNAFastaFile" ncrna.fa)${gene_biotype}.fa"

        echo "Output rRNA FASTA file path: $PathrRNAFastaFile"

        echo "Subsetting ncRNA FASTA file..."

        awk '
            BEGIN { keep=0 } 
            /^>/ { 
                keep = ($0 ~ /gene_biotype:rRNA/) ? 1 : 0 
            } 
            keep { print }
        ' "$PathncRNAFastaFile" > "$PathrRNAFastaFile"

        if [[ -f "$PathrRNAFastaFile" ]]; then
        log " - The $gene_biotype FASTA file have been created successfully"
            
        else
            log " - ERROR: $gene_biotype FASTA file not found."
            exit 1
        fi
    
        /usr/bin/time saw makeRef \
            --mode=STAR \
            --genome="${referenceDir}/makeRef" \
            --fasta="${genomeFastaFiles}" \
            --gtf="${annotationCheckGTFFiles}" \
            --rRNA-fasta ${PathrRNAFastaFile} \
            --threads-num="$threads"

        log " - Saw index complete."            

    fi

else

    log " - Creating Saw index..."
    
    /usr/bin/time saw makeRef \
        --mode=STAR \
        --genome="${referenceDir}/makeRef" \
        --fasta="${genomeFastaFiles}" \
        --gtf="${annotationCheckGTFFiles}" \
        --threads-num="$threads"

    log " - Saw index complete."


fi

rm -rf "${referenceDir}/genome/" "$referenceDir/genes/"
