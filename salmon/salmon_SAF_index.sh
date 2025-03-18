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
        echo "Cannot find conda initialization script."
        exit 1
    fi

    conda activate "$condaEnv" || { echo "Failed to activate conda environment: $condaEnv"; exit 1; }
}


deactivate_conda() {
    # Ensure the conda initialization script is sourced
    if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        . "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        . "/opt/conda/etc/profile.d/conda.sh"
    else
        echo "Cannot find conda initialization script."
        exit 1
    fi

    conda deactivate  || { echo "Failed to deactivate conda environment: $condaEnv"; exit 1; }
}

gtf2tx2gene() {
  # Check that an input file is provided
  if [ -z "$1" ]; then
    echo "Usage: gtf2tx2gene <file.gtf(.gz)> [output_file.tsv]"
    return 1
  fi

  local input="$1"
  local output="${2:-tx2gene_ensemble.tsv}"

  zless "$input" | grep -v "^#" | \
  awk 'BEGIN {
    FS="\t"; OFS="\t"
  }
  $3 == "transcript" {
    # Split the 9th column (attributes) by ";" to extract fields
    n = split($9, attr, ";");
    gene_id = "";
    transcript_id = "";
    
    # Loop through the attributes to extract gene_id and transcript_id
    for(i = 1; i <= n; i++){
      gsub(/^ +| +$/, "", attr[i]);
      
      # If the attribute starts with gene_id
      if(attr[i] ~ /^gene_id/){
        split(attr[i], tmp, " ");
        # tmp[2] contains the value (with quotes, which are removed)
        gene_id = tmp[2];
        gsub(/"/, "", gene_id);
      }
      # If the attribute starts with transcript_id
      else if(attr[i] ~ /^transcript_id/){
        split(attr[i], tmp, " ");
        transcript_id = tmp[2];
        gsub(/"/, "", transcript_id);
      }
    }

    # Print only if both fields were found
    if(gene_id != "" && transcript_id != ""){
      print transcript_id, gene_id
    }
  }' | sort -u > "$output"
  
  echo "tx2gene Files generated: $output"
}


usage() {
    echo "$(date) usage: $0 -PathOutputReference <path> -species <name> [-threads <num> -fasta <path> -gtf <path>]
    -PathOutputReference : Directory where reference will be created (required)
    -species : Organism type of sample, usually referring to species. (required)
    -condaEnv: CondaEnv with salmon (required)
    -threads : the number of threads to be used in running the pipeline (default: 4) (optional)
    -PathGenFastaFile : path to custom .fa file (optional)
    -PathTranFastaFile : path to custom .fa file (optional)
    -tx2gene : Create tx2gene file (optional)
    -PathGtfFile : path to custom .gtf file (optional)
    
    "
    exit 1
}
# Initialize variables
threads=4
LC_TIME=fr_FR.UTF-8
condaEnv=salmon

if [[ $# -lt 6 || $# -gt 12 ]]; then
    usage
fi

while [[ -n "$1" ]]; do
    case "$1" in
        -PathOutputReference) PathOutputReference="$2"; shift ;;
        -species) species="$2"; shift ;;
        -threads) threads="$2"; shift ;;
        -condaEnv) condaEnv="$2"; shift ;;
        -PathGenFastaFile) PathGenFastaFile="$2"; shift ;;
        -PathTranFastaFile) PathTranFastaFile="$2"; shift ;;
        -tx2gene) tx2gene_flag=1 ;;  
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
if [[ -z "$species" && (-z "$PathGenFastaFile" || -z "$PathTranFastaFile") ]]; then
    echo "$(date) - Either species name or both DNA/cDNA FASTA files must be provided."
    exit 1
fi

# If fasta and gtf are provided, validate them
if [[ -n "$PathGenFastaFile" && -n "$PathTranFastaFile" ]]; then
    if [[ ! -f "$PathGenFastaFile" ]]; then
        echo "$(date) - Provided DNA FASTA file does not exist: $PathGenFastaFile"
        exit 1
    fi
    if [[ ! -f "$PathTranFastaFile" ]]; then
        echo "$(date) - Provided cDNA FASTA file does not exist: $PathTranFastaFile"
        exit 1
    fi
else
    # If only one of the custom files is provided, exit with error
    if [[ -n "$PathGenFastaFile" || -n "$PathTranFastaFile" ]]; then
        echo "$(date) - Both fasta files must be provided together."
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
    echo "Directory $referenceDir already exist. Remove $referenceDir please or wait one minutes"
    exit 1
else
    echo "Directory $referenceDir does not exist. Creating directory."
    mkdir -p "$referenceDir"
    echo "Created new subdirectory: $referenceDir"
fi

# Use referenceDir for further operations
echo "Using $referenceDir for your operations."


# Create necessary directories
mkdir -p "$referenceDir/genome/cdna" "$referenceDir/genome/dna" "$referenceDir/genes"

cd "$referenceDir" || {
    echo "$(date) - Failed to change directory to $referenceDir"
    exit 1
}



# Logic for custom files
if [[ -n "$PathGenFastaFile" && -n "$PathTranFastaFile" ]]; then
    cp "$PathGenFastaFile" "$referenceDir/genome/dna"
    cp "$PathTranFastaFile" "$referenceDir/genome/cdna"
    
    if [[ "$PathGenFastaFile" == *.gz ]]; then
        
        genomeFastaFiles=$(find "${referenceDir}/genome/dna/" -type f -name "$(basename "$PathGenFastaFile" .gz)" | head -n 1)
    else
    gzip -f "$referenceDir/genome/dna/$(basename "$PathGenFastaFile")"
        genomeFastaFiles=$(find "${referenceDir}/genome/dna/" -type f -name "$(basename "$PathGenFastaFile")" | head -n 1)
    fi

    if [[ "$PathTranFastaFile" == *.gz ]]; then
  
        transcriptomeFastaFiles=$(find "${referenceDir}/genome/cdna/" -type f -name "$(basename "$PathTranFastaFile" .gz)" | head -n 1)
    else
        gzip -f "$referenceDir/genome/cdna/$(basename "$PathTranFastaFile")"
        transcriptomeFastaFiles=$(find "${referenceDir}/genome/cdna/" -type f -name "$(basename "$PathTranFastaFile")" | head -n 1)
    fi

    echo "$(date) - Creating decoys files ..."
    

    grep "^>" <(gunzip -c $genomeFastaFiles) | cut -d " " -f 1 > ${referenceDir}/genome/decoys.txt
    sed -i.bak -e 's/>//g' ${PathOutputReference}/decoys.txt


    cat ${transcriptomeFastaFiles} ${genomeFastaFiles} > ${referenceDir}/genome/gentrome.fa.gz


    echo "$(date) - Creating Salmon Index..."

    activate_conda

    /usr/bin/time salmon index -t ${referenceDir}/genome/gentrome.fa.gz -d ${referenceDir}/genome/decoys.txt -p $threads -i ${referenceDir}/salmon_index

    deactivate_conda

    rm -rf ${referenceDir}/genome/


else


    ensembl_species="${species,,}"

    request_dna_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/dna/"
    request_cdna_fasta_fd="ftp://ftp.ensembl.org/pub/current_fasta/${ensembl_species}/cdna/"
    genome_filename="*.dna.primary_assembly.fa.gz"
    transcriptome_filename="*.cdna.all.fa.gz"     


    echo "$(date) - Trying to download DNA primary assembly fasta file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genome/dna/" -A "$genome_filename" "$request_dna_fasta_fd" || {
        echo "$(date) - Primary assembly file not found. Trying to download DNA toplevel fasta file..."
        genome_filename="*.dna.toplevel.fa.gz"
        wget -r -np -nd -q -P "$referenceDir/genome/dna/" -A "$genome_filename" "$request_dna_fasta_fd" || {
            echo "$(date) - Failed to download DNA file from Ensembl."
            exit 1

        }
    }


    echo "$(date) - Trying to download cDNA fasta file from Ensembl..."
    wget -r -np -nd -q -P "$referenceDir/genome/cdna/" -A "${transcriptome_filename}" "$request_cdna_fasta_fd" || {
        echo "$(date) - Failed to download cDNA fasta file from Ensembl."
        exit 1
    }

    echo "$(date) - The DNA and cDNA fasta files have been downloaded successfully"


    genomeFastaFiles=$(find "${referenceDir}/genome/dna/" -type f -regex "*.dna.\(primary_assembly\|toplevel\).fa.gz" | head -n 1)
    transcriptomeFastaFiles=$(find "${referenceDir}/genome/cdna/" -type f -regex "*.cdna.all.fa.gz" | head -n 1)

    echo "$(date) - Creating decoys files ..."
    

    grep "^>" <(gunzip -c $genomeFastaFiles) | cut -d " " -f 1 > ${referenceDir}/genome/decoys.txt
    sed -i.bak -e 's/>//g' ${PathOutputReference}/decoys.txt


    cat ${transcriptomeFastaFiles} ${genomeFastaFiles} > ${referenceDir}/genome/gentrome.fa.gz


    echo "$(date) - Creating Salmon Index..."

    activate_conda

    /usr/bin/time salmon index -t ${referenceDir}/genome/gentrome.fa.gz -d ${referenceDir}/genome/decoys.txt -p $threads -i ${referenceDir}/salmon_index

    deactivate_conda

    rm -rf ${referenceDir}/genome/

fi


if [[ "$tx2gene_flag" -eq 1 ]]; then  # Run ONLY if -tx2gene flag is provided
    if [[ -f "$PathGtfFile" && ( "$PathGtfFile" == *.gtf || "$PathGtfFile" == *.gtf.gz ) ]]; then
        echo "$(date) - Creating tx2gene file from $PathGtfFile"
        mkdir -p "$referenceDir/genes/"
        cp "$PathGtfFile" "$referenceDir/genes/"
        # Compress if the file is a plain .gtf
        if [[ "$PathGtfFile" == *.gtf ]]; then
            gzip "${referenceDir}/genes/$(basename "$PathGtfFile")"
        fi
        # Find the GTF file (now guaranteed to be .gtf.gz)
        pathFileGTF=$(find "$referenceDir/genes/" -type f -name "*.gtf.gz" | head -n 1)
        if [[ -f "$pathFileGTF" ]]; then
            /usr/bin/time gtf2tx2gene "$pathFileGTF" ${referenceDir}/tx2gene.tsv
            echo "$(date) - ${referenceDir}/tx2gene.tsv file successfuly created "
            rm -rf ${referenceDir}/genes/
        else
            echo "$(date) - ERROR: GTF file not found after copying."
            exit 1
        fi
    else 
        echo "$(date) - Provided GTF file is invalid or missing."
        echo "$(date) - Downloading GTF file from Ensembl..."
        ensembl_species="${species,,}"  # Ensure lowercase species name
        request_gtf_fd="ftp://ftp.ensembl.org/pub/current_gtf/${ensembl_species}/"
        mkdir -p "$referenceDir/genes/"
        # Download the latest GTF
        wget -r -np -nd -q -P "$referenceDir/genes/" -A "*[0-9].gtf.gz" "$request_gtf_fd" || {
            echo "$(date) - ERROR: Failed to download GTF from Ensembl."
            exit 1
        }
        # Verify downloaded file
        pathFileGTF=$(find "$referenceDir/genes/" -type f -name "*.gtf.gz" | head -n 1)
        if [[ -f "$pathFileGTF" ]]; then
            /usr/bin/time gtf2tx2gene "$pathFileGTF" ${referenceDir}/tx2gene.tsv
            echo "$(date) - ${referenceDir}/tx2gene.tsv file successfuly created "
            rm -rf ${referenceDir}/genes/
        else
            echo "$(date) - ERROR: No GTF file found after download."
            exit 1
        fi
    fi
fi

    






