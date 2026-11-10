# Metagenomics Analysis Pipeline

This repository contains a collection of production-ready SLURM shell scripts for processing, filtering, assembling, and analyzing high-throughput metagenomic sequencing data on a high-performance computing (HPC) cluster.

---

## Workflow Overview

```text
Raw Reads 
  ├──> 1. Raw Quality Control (FastQC / MultiQC)
  ├──> 2. Trimming & Filtering (fastp)
  ├──> 3. Host Decontamination (Bowtie2 vs. GRCh38)
        ├──> 4a. De Novo Assembly (MEGAHIT)
        ├──> 4b. Taxonomic Profiling (Kraken2)
        └──> 4c. Reference Alignment (Minimap2)
```

---

## 1. Raw Quality Control (`qc.sh`)
Assesses the quality of raw paired-end sequencing reads prior to any downstream processing.

```bash
#!/bin/bash
#SBATCH --job-name=mg_qc
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/resources/datasets/metagenomics/raw_reads"
OUTPUT_DIR="/mnt/roche/${USER}_roche_project/results/qc"

module load fastqc multiqc

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "Job started on: $(date)"

# Run FastQC on all raw files
fastqc -o "$OUTPUT_DIR" "$INPUT_DIR"/*.fastq.gz

# Aggregate reports with MultiQC
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc_report"

echo "Job completed successfully on: $(date)"
```

---

## 2. Quality Trimming & Filtering (`fastp.sh`)
Removes sequencing adapters and low-quality bases using `fastp`.

```bash
#!/bin/bash
#SBATCH --job-name=mg_qc_fastp
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/resources/datasets/metagenomics/raw_reads"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
CLEAN_DIR="$BASE_OUT/cleaned_reads"
QC_DIR="$BASE_OUT/qc_reports"

module load fastp fastqc multiqc

mkdir -p "$CLEAN_DIR" "$QC_DIR" logs

for r1 in "$INPUT_DIR"/*_1.fastq.gz; do
    [ -e "$r1" ] || continue
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    base=$(basename "$r1" _1.fastq.gz)
    
    if [ -f "$r2" ]; then
        fastp \
            -i "$r1" -I "$r2" \
            -o "$CLEAN_DIR/${base}_1_clean.fastq.gz" \
            -O "$CLEAN_DIR/${base}_2_clean.fastq.gz" \
            -j "$QC_DIR/${base}_fastp.json" \
            -h "$QC_DIR/${base}_fastp.html" \
            --detect_adapter_for_pe
    fi
done

fastqc -o "$QC_DIR" "$CLEAN_DIR"/*.fastq.gz
multiqc "$QC_DIR" -o "$QC_DIR/multiqc_report"
```

---

## 3. Host Decontamination (`host_removal.sh`)
Aligns clean reads against the human reference genome (`GRCh38_noalt_as`) using `Bowtie2` and extracts unmapped (host-free) paired reads via `samtools`.

```bash
#!/bin/bash
#SBATCH --job-name=mg_host_removal
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/${USER}_roche_project/results/cleaned_reads"
HUMAN_INDEX="/mnt/roche/alfred_roche_project/databases/human/GRCh38_noalt_as"

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
HOSTFREE_DIR="$BASE_OUT/host_free_reads"
QC_DIR="$BASE_OUT/host_removal_qc"

module load bowtie2 samtools fastqc multiqc

mkdir -p "$HOSTFREE_DIR" "$QC_DIR" logs

for r1_clean in "$INPUT_DIR"/*_1_clean.fastq.gz; do
    [ -e "$r1_clean" ] || continue
    r2_clean="${r1_clean/_1_clean.fastq.gz/_2_clean.fastq.gz}"
    base=$(basename "$r1_clean" _1_clean.fastq.gz)
    
    if [ -f "$r2_clean" ]; then
        bowtie2 -x "$HUMAN_INDEX" \
            -1 "$r1_clean" -2 "$r2_clean" \
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

fastqc -t 4 -o "$QC_DIR" "$HOSTFREE_DIR"/*.fastq.gz
multiqc "$QC_DIR" -o "$QC_DIR/multiqc_report"
```

---

## 4. De Novo Assembly (`assembly.sh`)
Assembles host-free metagenomic reads into contigs using **MEGAHIT**.

```bash
#!/bin/bash
#SBATCH --job-name=mg_assembly
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=5:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
BASE_OUT="/mnt/roche/${USER}_roche_project/results"
ASSEMBLY_DIR="$BASE_OUT/assemblies"
MEGAHIT_BIN="/software/megahit/1.2.9/bin/megahit"

mkdir -p "$ASSEMBLY_DIR" logs

for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
    [ -e "$r1_hostfree" ] || continue
    r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
    base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
    
    if [ -f "$r2_hostfree" ]; then
        "$MEGAHIT_BIN" \
            -1 "$r1_hostfree" -2 "$r2_hostfree" \
            -t 16 \
            --out-dir "$ASSEMBLY_DIR/$base" \
            --out-prefix "$base"
    fi
done
```

---

## 5. Taxonomic Profiling (`kraken2.sh`)
Classifies host-free reads using **Kraken2** against a standard taxonomic database.

```bash
#!/bin/bash
#SBATCH --job-name=mg_kraken2
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
KRAKEN_DB="/mnt/roche/alfred_roche_project/databases/krakenDB"

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
KRAKEN_OUT="$BASE_OUT/kraken2_reports"

module load kraken2

mkdir -p "$KRAKEN_OUT" logs

for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
    [ -e "$r1_hostfree" ] || continue
    r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
    base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
    
    if [ -f "$r2_hostfree" ]; then
        kraken2 \
            --db "$KRAKEN_DB" \
            --threads 4 \
            --paired \
            --output "$KRAKEN_OUT/${base}.kraken" \
            --report "$KRAKEN_OUT/${base}.report.txt" \
            "$r1_hostfree" "$r2_hostfree"
    fi
done
```

---

## 6. Reference-Based Alignment (`mapping.sh`)
Maps host-free reads against custom reference genomes using **minimap2** and sorts/indexes the resulting BAM files.

```bash
#!/bin/bash
#SBATCH --job-name=mg_minimap2
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

INPUT_DIR="/mnt/roche/${USER}_roche_project/results/host_free_reads"
REF_DIR="/mnt/roche/resources/metagenomics/databases/ref"

ACCESSION_ID="YOUR_ACCESSION_ID_HERE"
REF_GENOME="$REF_DIR/${ACCESSION_ID}.fa"

BASE_OUT="/mnt/roche/${USER}_roche_project/results"
MAP_DIR="$BASE_OUT/mapping_${ACCESSION_ID}"
MINIMAP_BIN="/software/minimap2/2.26/minimap2"

module load samtools

mkdir -p "$MAP_DIR" logs

for r1_hostfree in "$INPUT_DIR"/*_1_hostfree.fastq.gz; do
    [ -e "$r1_hostfree" ] || continue
    r2_hostfree="${r1_hostfree/_1_hostfree.fastq.gz/_2_hostfree.fastq.gz}"
    base=$(basename "$r1_hostfree" _1_hostfree.fastq.gz)
    
    if [ -f "$r2_hostfree" ]; then
        "$MINIMAP_BIN" -ax sr -t 8 "$REF_GENOME" "$r1_hostfree" "$r2_hostfree" | \
        samtools view -b -@ 2 - | \
        samtools sort -@ 4 -o "$MAP_DIR/${base}_vs_${ACCESSION_ID}.sorted.bam"
        
        samtools index "$MAP_DIR/${base}_vs_${ACCESSION_ID}.sorted.bam"
    fi
done
```
