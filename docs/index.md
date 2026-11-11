# Metagenomics Analysis Pipeline

This repository contains a collection of production-ready SLURM shell scripts for processing, filtering, assembling, and analyzing high-throughput metagenomic sequencing data on a high-performance computing (HPC) cluster.

## Getting Started & Workspace Guidelines

To maintain an organized environment and avoid permission conflicts:

1. **Personal Workspace:** Every user must work strictly inside their individual workspace directory:
   

   ```
   /mnt/roche/${USER}_roche_project
   ```

2. **Shared Resources:** Raw sample datasets and pre-downloaded reference databases are already provided and maintained in the central shared folder:
   

   ```
   /mnt/roche/resources/datasets/metagenomics
   ```

## Expected Project File Structure

Inside your personal workspace (`/mnt/roche/${USER}_roche_project`), your results and project layout should follow this structure:

```
/mnt/roche/${USER}_roche_project/
├── scripts/
│   ├── qc.sh                  # 1. Raw Quality Control
│   ├── trimming.sh            # 2. Trimming & Filtering
│   ├── host_removal.sh        # 3. Host Decontamination
│   ├── kraken2.sh             # 4a. Taxonomic Profiling
│   ├── run_krona.sh           # 4b. Taxonomic Visualization
│   ├── mapping.sh             # 5. Reference-Based Alignment
│   ├── denovo.sh              # 6a. De Novo Assembly
│   └── quast.sh               # 6b. Assembly Evaluation
│
├── results/
│   ├── qc/                       # Raw FastQC & MultiQC outputs
│   ├── cleaned_reads/            # Fastp trimmed and filtered reads
│   ├── qc_reports/               # Fastp JSON/HTML and MultiQC reports
│   ├── host_free_reads/          # Host-decontaminated paired FASTQ files
│   ├── host_removal_qc/          # FastQC & MultiQC on host-free reads
│   ├── kraken2_reports/          # Kraken2 classification reports (.report.txt)
│   ├── krona_plots/              # Interactive Krona HTML visualizations
│   ├── mapping_<ACCESSION_ID>/   # Sorted BAM files and indices
│   ├── assemblies/               # MEGAHIT contig directories per sample
│   └── quast_reports/            # MetaQUAST evaluation reports
│
└── logs/                         # SLURM job standard output and error logs

```

---

## Initializing the Workspace & Running the Pipeline
Before running any pipeline steps, set up your workspace directory structure and populate your scripts/ folder from the shared resource directory:

```{bash}
# 1. Navigate to your project root workspace
cd /mnt/roche/${USER}_roche_project

# 2. Create your project directories in your workspace
mkdir -p /mnt/roche/${USER}_roche_project/{scripts,results,logs}

# 3. Copy the pre-prepared SLURM scripts into your scripts folder
cp /mnt/roche/resources/datasets/metagenomics/scripts/*.sh /mnt/roche/${USER}_roche_project/scripts/

# 4. Execute pipeline steps directly from the root workspace using SLURM
sbatch scripts/qc.sh
sbatch scripts/trimming.sh
sbatch scripts/host_removal.sh
sbatch scripts/kraken2.sh
sbatch scripts/run_krona.sh
sbatch scripts/mapping.sh
sbatch scripts/denovo.sh
sbatch scripts/quast.sh

# 5. Monitor your job queue status
squeue -u $USER
```

---

## Raw Quality Control (`qc.sh`)

We perform initial raw quality control to assess the structural integrity of the sequencing data, identify adapter contamination, detect per-base sequence quality degradation, and flag any major sequencing artifacts or GC-bias before spending computational resources on downstream analyses.

```bash
#!/bin/bash
#SBATCH --job-name=mg_qc
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%j_qc.out
#SBATCH --error=logs/%j_qc.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/resources/datasets/metagenomics/raw_reads"
OUTPUT_DIR="/mnt/roche/${USER}_roche_project/results/qc"

# ===================== ENVIRONMENT SETUP =====================
module load fastqc multiqc


# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# ===================== RUN FASTQC =====================
echo "Running FastQC..."

# Adjust the file extension wildcard if your files end in .fq.gz or .fastq
fastqc -o "$OUTPUT_DIR" "$INPUT_DIR"/*.fastq.gz

echo "FastQC finished."

# ===================== RUN MULTIQC =====================
echo "Running MultiQC to aggregate reports..."

# Multiqc will scan the output directory for fastqc zip/html files
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc_report"

echo "Job completed successfully on: $(date)"
```

After running MultiQC, check per-base sequence quality scores, duplication levels, and adapter-content metrics across samples to determine the stringency needed for subsequent trimming and filtering.

---

## Quality Trimming & Filtering (`trimming.sh`)

Trimming and filtering are essential to remove low-quality bases, ambiguous nucleotides, and remaining adapter sequences that can otherwise distort downstream sequence alignments, inflate assembly fragmentation, or introduce false variants.

```bash
#!/bin/bash
#SBATCH --job-name=mg_qc_fastp
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j_trimming.out
#SBATCH --error=logs/%j_trimming.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/resources/datasets/metagenomics/raw_reads"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
CLEAN_DIR="$BASE_OUT/cleaned_reads"
QC_DIR="$BASE_OUT/qc_reports"

# ===================== ENVIRONMENT SETUP =====================
module load fastp fastqc multiqc

# Create necessary directories
mkdir -p "$CLEAN_DIR"
mkdir -p "$QC_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory: $INPUT_DIR"

# ===================== RUN FASTP =====================
echo "Running fastp quality filtering and adapter trimming..."

# Loop through paired-end reads (adjust suffix if your files use _R1.fastq.gz or .fq.gz)
for r1 in "$INPUT_DIR"/*_1.fastq.gz; do
    # Check if the file exists to prevent glob literal errors if no match
    [ -e "$r1" ] || continue
    
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    base=$(basename "$r1" _1.fastq.gz)
    
    if [ -f "$r2" ]; then
        echo "Processing sample: $base"
        fastp \
            -i "$r1" \
            -I "$r2" \
            -o "$CLEAN_DIR/${base}_1_clean.fastq.gz" \
            -O "$CLEAN_DIR/${base}_2_clean.fastq.gz" \
            -j "$QC_DIR/${base}_fastp.json" \
            -h "$QC_DIR/${base}_fastp.html" \
            --detect_adapter_for_pe
    fi
done

echo "fastp finished."

# ===================== RUN FASTQC =====================
echo "Running FastQC on cleaned reads..."
fastqc -o "$QC_DIR" "$CLEAN_DIR"/*.fastq.gz

echo "FastQC finished."

# ===================== RUN MULTIQC =====================
echo "Running MultiQC to aggregate fastp and FastQC reports..."
multiqc "$QC_DIR" -o "$QC_DIR/multiqc_report"

echo "Job completed successfully on: $(date)"
```

Review the fastp JSON summary and the post-filtering FastQC reports to verify that low-quality tails have been successfully removed while ensuring an acceptable proportion of clean, usable read pairs remains for downstream applications.

---

## Host Decontamination (`host_removal.sh`)

In clinical or host-associated metagenomics, a significant portion of raw reads often originate from the host organism rather than the target microbial or viral community; removing these host sequences reduces computational overhead and prevents misattribution during profiling. This step aligns clean reads against the human reference genome (`GRCh38_noalt_as`) using `Bowtie2` and extracts unmapped (host-free) paired reads via `samtools`.

```bash
#!/bin/bash
#SBATCH --job-name=mg_host_removal
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j_host_removal.out
#SBATCH --error=logs/%j_host_removal.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/${USER}_roche_project/results/cleaned_reads"
HUMAN_INDEX="/mnt/roche/${USER}_roche_project/databases/human/GRCh38_noalt_as"
BOWTIE2_BIN="/software/bowtie2/2.4.2/bowtie2"

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
HOSTFREE_DIR="$BASE_OUT/host_free_reads"
QC_DIR="$BASE_OUT/host_removal_qc"

# ===================== ENVIRONMENT SETUP =====================
module load samtools fastqc multiqc

# Create necessary directories
mkdir -p "$HOSTFREE_DIR"
mkdir -p "$QC_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory (from fastp): $INPUT_DIR"
echo "Human reference index: $HUMAN_INDEX"
echo "Output host-free directory: $HOSTFREE_DIR"

# ===================== RUN BOWTIE2 HOST REMOVAL =====================
echo "Running Bowtie2 to filter out human host reads..."

for r1_clean in "$INPUT_DIR"/*_1_clean.fastq.gz; do
    # Check if files exist to prevent glob errors
    [ -e "$r1_clean" ] || continue
    
    r2_clean="${r1_clean/_1_clean.fastq.gz/_2_clean.fastq.gz}"
    base=$(basename "$r1_clean" _1_clean.fastq.gz)
    
    if [ -f "$r2_clean" ]; then
        echo "Processing sample: $base"
        
        # 1. Align reads against human genome index with Bowtie2 using direct path
        # 2. Use samtools view -f 12 to filter and keep only unmapped read pairs (both mates unmapped)
        # 3. Sort by name and convert back to paired FASTQ files for downstream metagenomic tools
        "$BOWTIE2_BIN" -x "$HUMAN_INDEX" \
            -1 "$r1_clean" \
            -2 "$r2_clean" \
            --threads 4 \
            2> "$QC_DIR/${base}_bowtie2.log" | \
        samtools view -b -f 12 - | \
        samtools sort -n -@ 4 - | \
        samtools fastq \
            -1 "$HOSTFREE_DIR/${base}_1_hostfree.fastq.gz" \
            -2 "$HOSTFREE_DIR/${base}_2_hostfree.fastq.gz" \
            -0 /dev/null -s /dev/null -
    fi
done

echo "Bowtie2 host removal finished."

# ===================== RUN FASTQC & MULTIQC =====================
echo "Running FastQC on host-free reads..."
fastqc -t 4 -o "$QC_DIR" "$HOSTFREE_DIR"/*.fastq.gz

echo "Running MultiQC to summarize host removal metrics..."
multiqc "$QC_DIR" -o "$QC_DIR/multiqc_report"

echo "Job completed successfully on: $(date)"
```

Inspect the Bowtie2 alignment logs and MultiQC summary metrics to determine the exact percentage of host reads filtered out per sample, ensuring that a sufficient non-host fraction remains for assembly and profiling.

---

## Taxonomic Profiling (`kraken2.sh`)

Kraken2 allows for rapid, k-mer-based taxonomic classification of sequencing reads against a reference database, giving us an immediate quantitative overview of the microbial composition, relative abundances, and presence of specific pathogens within the sample.

```bash
#!/bin/bash
#SBATCH --job-name=mg_kraken2
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j_kraken.out
#SBATCH --error=logs/%j_kraken.err

# ===================== EDIT HERE =====================
# Input points to the host-free reads generated from the previous step
INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
KRAKEN_DB="/mnt/roche/resources/datasets/metagenomics/databases/krakenDB"

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
KRAKEN_OUT="$BASE_OUT/kraken2_reports"

# ===================== ENVIRONMENT SETUP =====================
module load kraken2

# Create necessary directories
mkdir -p "$KRAKEN_OUT"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory (host-free reads): $INPUT_DIR"
echo "Kraken2 database: $KRAKEN_DB"
echo "Output directory: $KRAKEN_OUT"

# ===================== RUN KRAKEN2 =====================
echo "Running Kraken2 taxonomic classification..."

for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
    # Check if files exist to prevent glob errors
    [ -e "$r1_hostfree" ] || continue
    
    r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
    base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
    
    if [ -f "$r2_hostfree" ]; then
        echo "Profiling sample with Kraken2: $base"
        
        # Run Kraken2 in paired-end mode
        # Note: If your database is extremely large and you encounter RAM limits, 
        # you can add --use-mmap to the command to map the database into memory.
        kraken2 \
            --db "$KRAKEN_DB" \
            --threads 4 \
            --paired \
            --output "$KRAKEN_OUT/${base}.kraken" \
            --report "$KRAKEN_OUT/${base}.report.txt" \
            "$r1_hostfree" "$r2_hostfree"
    fi
done

echo "Kraken2 classification finished."
```

Examine the Kraken2 report files to check the classification rates, verify the proportion of unclassified versus classified reads, and identify dominant taxa or specific pathogen signatures present in each sample.


### Visualise Kraken reports

Visualizing taxonomic profiles through hierarchical interactive charts (Krona) allows researchers to intuitively explore microbial diversity across multiple taxonomic ranks—from phylum down to species—making it easy to spot dominant populations, rare taxa, or unexpected shifts in community composition at a glance.

```bash
#!/bin/bash
#SBATCH --job-name=mg_krona
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%j_krona.out
#SBATCH --error=logs/%j_krona.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/${USER}_roche_project/results/kraken2_reports"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
KRONA_DIR="$BASE_OUT/krona_plots"
KRAKENTOOLS_DIR="/software/krakentools/1.2.1"

# ===================== ENVIRONMENT SETUP =====================
module load krona

# Create necessary directories
mkdir -p "$KRONA_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory (Kraken reports): $INPUT_DIR"
echo "Output directory (Krona plots): $KRONA_DIR"

# ===================== RUN KRONATOOLS VISUALIZATION =====================
echo "Generating Krona interactive charts from Kraken2 reports..."

for report in "$INPUT_DIR"/*.report.txt; do
    # Check if files exist to prevent glob errors
    [ -e "$report" ] || continue
    base=$(basename "$report" .report.txt)
    
    echo "Processing sample: $base"
    
    # Step 1: Convert Kraken report to Krona format using KrakenTools
    python3 "$KRAKENTOOLS_DIR/kreport2krona.py" \
        --report "$report" \
        --output "$KRONA_DIR/${base}.krona.txt"
        
    # Step 2: Generate interactive HTML plot using KronaTools
    ktImportText "$KRONA_DIR/${base}.krona.txt" \
        -o "$KRONA_DIR/${base}_krona.html"
        
    # Optional: Remove intermediate text file
    rm "$KRONA_DIR/${base}.krona.txt"
done

echo "Krona visualization finished successfully on: $(date)"
```

Open the generated `.html` files in any web browser to interactively drill down through taxonomic branches, evaluate the relative proportions of major groups, check for unexpected high-abundance species that could point to reagent or sample contamination, and assess the overall depth of taxonomic resolution achieved.

---

## Reference-Based Alignment (`mapping.sh`)
Reference-based alignment maps our host-free reads directly against a specific reference genome of interest, enabling deep variant calling, consensus generation, coverage depth assessment, and micro-diversity tracking for target organisms.

```bash
#!/bin/bash
#SBATCH --job-name=mg_minimap2
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%j_mapping.out
#SBATCH --error=logs/%j_mapping.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
REF_DIR="/mnt/roche/resources/datasets/metagenomics/databases/refs"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
MAP_BASE_DIR="$BASE_OUT/reference_mapping"

# ===================== ENVIRONMENT SETUP =====================
module load minimap2 samtools

# Create necessary directories
mkdir -p "$MAP_BASE_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input reads directory: $INPUT_DIR"
echo "Reference directory: $REF_DIR"
echo "Output mapping directory: $MAP_BASE_DIR"

# ===================== RUN MINIMAP2 ACROSS REFERENCES =====================
echo "Scanning reference directory for target genomes..."

for ref_file in "$REF_DIR"/*.fasta; do
    # Check if reference files exist to prevent glob errors
    [ -e "$ref_file" ] || continue
    
    # Extract accession ID from filename (e.g., NC_001348 from NC_001348.fasta)
    accession=$(basename "$ref_file" .fasta)
    map_dir="$MAP_BASE_DIR/${accession}"
    mkdir -p "$map_dir"
    
    echo "=================================================="
    echo "Target Reference: $accession ($ref_file)"
    echo "=================================================="
    
    for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
        [ -e "$r1_hostfree" ] || continue
        
        r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
        base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
        
        if [ -f "$r2_hostfree" ]; then
            echo "Mapping sample $base against $accession..."
            
            # Minimap2 command:
            # -ax sr: Preset for accurate short genomic paired-end reads
            # -t 8: Matches --cpus-per-task=8
            # Samtools: Converts SAM to sorted BAM and indexes it
            minimap2 -ax sr -t 8 "$ref_file" "$r1_hostfree" "$r2_hostfree" | \
            samtools view -b -@ 2 - | \
            samtools sort -@ 4 -o "$map_dir/${base}_vs_${accession}.sorted.bam"
            
            # Index the sorted BAM file
            samtools index "$map_dir/${base}_vs_${accession}.sorted.bam"
        fi
    done
done

echo "All minimap2 reference mapping and sorting tasks completed successfully on: $(date)"
```

Inspect the generated sorted BAM files and alignment statistics to check genome-wide coverage breadth, mean sequencing depth, and mapping rates against the target reference accession.

## De Novo Assembly (`assembly.sh`)

De novo assembly reconstructs overlapping reads into longer contiguous sequences (contigs), allowing us to recover novel viral genomes, bacterial strains, or functional genes without relying strictly on existing reference databases.
Assembles host-free metagenomic reads into contigs using **MEGAHIT**.

```bash
#!/bin/bash
#SBATCH --job-name=mg_assembly
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=logs/%j_denovo_assembly.out
#SBATCH --error=logs/%j_denovo_assembly.err

# ===================== EDIT HERE =====================
INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
ASSEMBLY_DIR="$BASE_OUT/assemblies"
MEGAHIT_BIN="/software/megahit/1.2.9/bin/megahit"

# ===================== ENVIRONMENT SETUP =====================
# Create necessary directories
mkdir -p "$ASSEMBLY_DIR"
mkdir -p logs

echo "Job started on: $(date)"
echo "Input directory (host-free reads): $INPUT_DIR"
echo "Output assembly directory: $ASSEMBLY_DIR"

# ===================== RUN MEGAHIT =====================
echo "Running MEGAHIT metagenomic assembly..."

for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
    # Check if files exist to prevent glob errors
    [ -e "$r1_hostfree" ] || continue
    
    r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
    base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
    
    if [ -f "$r2_hostfree" ]; then
        echo "Assembling sample: $base"
        
        # MEGAHIT execution using direct path
        "$MEGAHIT_BIN" \
            -1 "$r1_hostfree" \
            -2 "$r2_hostfree" \
            -t 16 \
            --out-dir "$ASSEMBLY_DIR/$base" \
            --out-prefix "$base"
    fi
done

echo "MEGAHIT assembly finished successfully on: $(date)"
```

Evaluate the assembly outputs by checking metrics like contig count, maximum length, and N50 to gauge the overall continuity and structural quality of the reconstructed genomic drafts.

---

### De Novo Assembly Evaluation (quast.sh)

Evaluating de novo assemblies using QUAST or MetaQUAST is essential because raw assembly outputs often contain fragmentation, chimeric contigs, or structural errors; these tools compute standardized metrics like N50, largest contig size, and total assembly length to objectively measure assembly contiguity and reliability.

```bash
#!/bin/bash
#SBATCH --job-name=mg_quast
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%j_quast.out
#SBATCH --error=logs/%j_quast.err

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
ASSEMBLY_DIR="$BASE_OUT/assemblies"
QUAST_OUT="$BASE_OUT/quast_reports"

module load quast

mkdir -p "$QUAST_OUT" logs

echo "Job started on: $(date)"

# Gather all contig files from individual assembly folders
contig_files=("$ASSEMBLY_DIR"/*/*.contigs.fa)

if [ ${#contig_files[@]} -gt 0 ]; then
    metaquast.py "${contig_files[@]}" \
        -o "$QUAST_OUT" \
        --threads 8
    echo "MetaQUAST evaluation finished."
else
    echo "No contig files found for QUAST evaluation."
fi

echo "Job completed successfully on: $(date)"
```

Examine the generated QUAST/MetaQUAST reports to check whether the N50 is sufficiently high and whether the total assembled length aligns with expectations, ensuring the assembly is robust enough for downstream binning, gene prediction, or functional annotation.

