#!/bin/bash
set -e
# Function to replace spaces with underscores and capitalize the first letter of each word
transform_and_capitalize() {
    echo "${1// /_}" | awk '{for(i=1;i<=NF;i++){$i=toupper(substr($i,1,1)) tolower(substr($i,2));}print}'
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

    conda activate "$condaEnv" || { log " - Failed to activate conda environment: $condaEnv"; exit 1; }
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

gtf2tx2gene() {
  # Check that an input file is provided
  if [ -z "$1" ]; then
    log " - Usage: gtf2tx2gene <file.gtf(.gz)> [output_file.tsv]"
    return 1
  fi

  local input="$1"
  local output="${2:-tx2gene_ensemble.tsv}"

  zless "$input" | grep -v "^#" | \
    awk 'BEGIN {
        FS="\t"; OFS="\t"
    }
    $3 == "transcript" {
        n = split($9, attr, ";");
        gene_id = ""; transcript_id = "";
        for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", attr[i]);
        if (attr[i] ~ /^gene_id/) {
            split(attr[i], tmp, " ");
            gene_id = tmp[2]; gsub(/"/, "", gene_id);
        } else if (attr[i] ~ /^transcript_id/) {
            split(attr[i], tmp, " ");
            transcript_id = tmp[2]; gsub(/"/, "", transcript_id);
        }
        }
        if (gene_id != "" && transcript_id != "" && !seen[transcript_id]++) {
        if (transcript_id ~ /^ENS/) {
            ens[transcript_id] = gene_id
        } else {
            other[transcript_id] = gene_id
        }
        }
    }
    END {
        for (t in ens)   print t, ens[t]
        for (t in other) print t, other[t]
    }' > "$output"
  
  echo "tx2gene Files generated: $output"
}
# Logging function for standardized message format
log() {
    echo "$(date) - $1"
}




usage() {
    log " usage: $0 -PathOutputReference <path> -speciesNames <name> -condaEnv <env> [-threads <num> -PathGenFastaFile <path> -PathTranFastaFile <path> -tx2gene -PathGtfFile <path>]
    
    Required arguments:
    -PathOutputReference : Directory where reference will be created.
    -speciesNames : Organism type of sample, usually referring to species.
    -condaEnv : Conda environment with Salmon installed.

    Optional arguments:
    -threads : Number of threads to use.
    -PathGenFastaFile : Path to genome FASTA file.
    -PathTranFastaFile : Path to transcript FASTA file.
    -tx2gene : Flag to create a tx2gene file.
    -PathGtfFile : Path to GTF annotation file to create tx2gene file.

    "
    exit 1
}
# Initialize variables
LC_TIME=fr_FR.UTF-8
condaEnv=salmon


while [[ -n "$1" ]]; do
    case "$1" in
        -PathOutputReference) PathOutputReference="$2"; shift ;;
        -speciesNames) speciesNames="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -condaEnv) condaEnv="$2"; shift ;;
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathTranFastaFile) PathTranFastaFile="$2"; shift ;;
        -tx2gene) tx2gene_flag=1 ;;  
        -PathGtfFile) PathGtfFile="$2"; shift ;;
        *) log " - Unknown option: $1"; usage ;;
    esac
    shift
done

# Check for mandatory arguments
if [[ -z "$PathOutputReference" || -z "$speciesNames" || -z "$condaEnv" ]]; then
    log " - Error: Missing required arguments."
    usage
fi


# If fasta and gtf are provided, validate them
if [[ -n "$PathGenFastaFile" && -n "$PathTranFastaFile" ]]; then
    if [[ ! -f "$PathGenFastaFile" ]]; then
        log " - Provided DNA FASTA file does not exist: $PathGenFastaFile"
        exit 1
    fi
    if [[ ! -f "$PathTranFastaFile" ]]; then
        log " - Provided cDNA FASTA file does not exist: $PathTranFastaFile"
        exit 1
    fi
else
    # If only one of the custom files is provided, exit with error
    if [[ -n "$PathGenFastaFile" || -n "$PathTranFastaFile" ]]; then
        log " - Both fasta files must be provided together."
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
referenceDir="${PathOutputReference}/${speciesNames}/${dateSuffix}"

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
mkdir -p "$referenceDir/genome/cdna" "$referenceDir/genome/dna" 

cd "$referenceDir" || {
    log " - Failed to change directory to $referenceDir"
    exit 1
}



# Logic for custom files
if [[ -n "$PathGenFastaFile" && -n "$PathTranFastaFile" ]]; then
    cp "$PathGenFastaFile" "$referenceDir/genome/dna"
    cp "$PathTranFastaFile" "$referenceDir/genome/cdna"
    
    if [[ "$PathGenFastaFile" != *.gz ]]; then
        gzip -f "$referenceDir/genome/dna/$(basename "$PathGenFastaFile")"   
    fi

    if [[ "$PathTranFastaFile" != *.gz ]]; then
        gzip -f "$referenceDir/genome/cdna/$(basename "$PathTranFastaFile")"
    fi

 

else


    ensembl_species="${speciesNames,,}"
    request_dna_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/dna/"
    request_cdna_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/cdna/"
    
        


    log " - Trying to download DNA primary assembly fasta file from Ensembl..."
    genome_filename="*.dna.primary_assembly.fa.gz"
    wget -r -np -nd -q -P "$referenceDir/genome/dna/" -A "$genome_filename" "$request_dna_fasta_fd" || {
        log " - Primary assembly file not found. Trying to download DNA toplevel fasta file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna/" -A "$genome_filename" "$request_dna_fasta_fd" || {
            log " - Failed to download DNA file from Ensembl."
            exit 1

        }
    }


    log " - Trying to download cDNA fasta file from Ensembl..."
    transcriptome_filename="*.cdna.all.fa.gz" 
    wget -r -np -nd -q -P "$referenceDir/genome/cdna/" -A "${transcriptome_filename}" "$request_cdna_fasta_fd" || {
        log " - Failed to download cDNA fasta file from Ensembl."
        exit 1
    }

    log " - The DNA and cDNA fasta files have been downloaded successfully"

fi

genomeFastaFiles=($(find "${referenceDir}/genome/dna/" -type f -name "*.fa.gz"))
transcriptomeFastaFiles=($(find "${referenceDir}/genome/cdna/" -type f -name "*.fa.gz"))

log " - Creating decoys files ..."

grep "^>" <(gunzip -c $genomeFastaFiles) | cut -d " " -f 1 > ${referenceDir}/genome/decoys.txt
sed -i.bak -e 's/>//g' ${referenceDir}/genome/decoys.txt

cat ${transcriptomeFastaFiles} ${genomeFastaFiles} > ${referenceDir}/genome/gentrome.fa.gz

log " - Creating Salmon index..."

activate_conda

/usr/bin/time salmon index \
  -t ${referenceDir}/genome/gentrome.fa.gz \
  -d ${referenceDir}/genome/decoys.txt \
  -p $threads \
  -i ${referenceDir}/makeRef

deactivate_conda

log " - Salmon index complete."
rm -rf ${referenceDir}/genome/

if [[ "$tx2gene_flag" -eq 1 ]]; then  # Run ONLY if -tx2gene flag is provided

    mkdir -p "$referenceDir/genes/"
    if [[ -f "$PathGtfFile" && ( "$PathGtfFile" == *.gtf || "$PathGtfFile" == *.gtf.gz ) ]]; then

        log " - Creating tx2gene file from $PathGtfFile"
        cp "$PathGtfFile" "$referenceDir/genes/"
        # Compress if the file is a plain .gtf
        if [[ "$PathGtfFile" == *.gtf ]]; then
            gzip "${referenceDir}/genes/$(basename "$PathGtfFile")"
        fi
        # Find the GTF file (now guaranteed to be .gtf.gz)
        pathFileGTF=$(find "$referenceDir/genes/" -type f -name "*.gtf.gz" | head -n 1)

        if [[ -f "$pathFileGTF" ]]; then

            gtf2tx2gene "$pathFileGTF" ${referenceDir}/tx2gene.tsv
            
        else
            log " - ERROR: GTF file not found after copying."
            exit 1
        fi
    else 
        log " - Provided GTF file is invalid or missing."
        log " - Downloading GTF file from Ensembl..."
        ensembl_species="${speciesNames,,}"  # Ensure lowercase species name
        request_gtf_fd="ftp://ftp.ensembl.org/pub/current_gtf/${ensembl_species}/"
        
        # Download the latest GTF
        wget -r -np -nd -q -P "$referenceDir/genes/" -A "*[0-9].gtf.gz" "$request_gtf_fd" || {
            log " - ERROR: Failed to download GTF from Ensembl."
            exit 1
        }
        log " - The gtf files have been downloaded successfully"
        # Verify downloaded file
        pathFileGTF=$(find "$referenceDir/genes/" -type f -regex ".*\.gtf\.gz" | head -n 1)

        if [[ -f "$pathFileGTF" ]]; then
            gtf2tx2gene "$pathFileGTF" ${referenceDir}/tx2gene.tsv
        else
            log " - ERROR: No GTF file found after download."
            exit 1
        fi
    fi
    log " - gtf2tx2gene complete."
    rm -rf ${referenceDir}/genes/
fi

    






