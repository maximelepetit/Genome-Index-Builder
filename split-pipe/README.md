# Saw 

---



## **1. Basic Usage (Ensembl Auto-Download)**
**Command**:
```bash

```
**What Happens**:

- Downloads genome (dna.primary_assembly.fa.gz) and transcriptome (cdna.all.fa.gz) from Ensembl.
- Generates Salmon index in /path/to/output/Gallus_gallus/YYYY_MM_DD/salmon_index.
- Does NOT generate tx2gene.tsv (no -tx2gene flag).
- 
**Validation**:
```bash
# Verify Salmon index files
ls -l /path/to/output/Gallus_gallus/YYYY_MM_DD/salmon_index
```

## **2. Basic Usage + tx2gene (Ensembl Auto-Download)**
**Command**:
```bash
./salmon_SAF_index.sh \
  -PathOutputReference "/path/to/output" \
  -species "Homo_sapiens" \
  -condaEnv salmon \
  -threads 16 \
  -tx2gene
```
**What Happens**:

- Downloads genome, transcriptome, and GTF from Ensembl.
- Generates Salmon index.
- Creates tx2gene.tsv from the downloaded GTF.

**Validation**:
```bash
# Check files
ls -l /path/to/output/Homo_sapiens/YYYY_MM_DD/
# Expected: salmon_index/, tx2gene.tsv
```

## **3. Custom FASTA Files **
**Command**:
```bash
./salmon_SAF_index.sh \
  -PathOutputReference "/path/to/output" \
  -species "Custom_species" \
  -condaEnv salmon \
  -threads 16 \
  -PathGenFastaFile "/path/to/genome.fa.gz" \
  -PathTranFastaFile "/path/to/transcriptome.fa.gz"
```
**What Happens**:

- Uses provided genome/transcriptome FASTA files.
- Generates Salmon index.
- Does NOT generate tx2gene.tsv.

**Validation**:
```bash
# Confirm index creation
ls -l /path/to/output/Custom_species/YYYY_MM_DD/salmon_index
```

## **4. Custom FASTA + Custom GTF + tx2gene**
**Command**:
```bash
./salmon_SAF_index.sh \
  -PathOutputReference "/path/to/output" \
  -species "Custom_species" \
  -condaEnv salmon \
  -threads 16 \
  -PathGenFastaFile "/path/to/genome.fa" \
  -PathTranFastaFile "/path/to/transcriptome.fa" \
  -tx2gene \
  -PathGtfFile "/path/to/annotations.gtf"
```
**What Happens**:

- Uses custom FASTA files (compresses them if needed).
- Creates tx2gene.tsv from the downloaded GTF.

**Validation**:
```bash
ls -l /path/to/output/Homo_sapiens/YYYY_MM_DD/
# Expected: salmon_index/, tx2gene.tsv,
```


## **Troubleshooting**
**Command**:
```bash
./salmon_SAF_index.sh \
  -PathOutputReference "/path/to/output" \
  -species "Invalid_species_123" \
  -condaEnv salmon
```
**Expected**: Fails with "Species name must contain only letters/underscores".

**Command**:
```bash
./salmon_SAF_index.sh \
  -PathOutputReference "/path/to/output" \
  -species "Mouse" \
  -PathGenFastaFile "/path/to/genome.fa"
```
**Expected**: Fails with "Both fasta files must be provided together".
