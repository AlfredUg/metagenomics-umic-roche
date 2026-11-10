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
├── scripts/                # Store your SLURM submission scripts here
├── logs/                   # Job output (.out) and error (.err) log files
└── results/
    ├── qc/                 # Raw FastQC and MultiQC reports
    ├── cleaned_reads/      # Trimming and filtering outputs (fastp)
    ├── host_free_reads/    # Host-decontaminated paired-end reads (Bowtie2)
    ├── assemblies/         # De novo contigs (MEGAHIT)
    ├── kraken2_reports/    # Taxonomic classification outputs
    └── mapping_*/          # Reference-based alignment BAM files (Minimap2)

```


---

## 1. Raw Quality Control (`qc.sh`)

We perform initial raw quality control to assess the structural integrity of the sequencing data, identify adapter contamination, detect per-base sequence quality degradation, and flag any major sequencing artifacts or GC-bias before spending computational resources on downstream analyses.

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

After running MultiQC, check per-base sequence quality scores, duplication levels, and adapter-content metrics across samples to determine the stringency needed for subsequent trimming and filtering.

---

## 2. Quality Trimming & Filtering (`fastp.sh`)

Trimming and filtering are essential to remove low-quality bases, ambiguous nucleotides, and remaining adapter sequences that can otherwise distort downstream sequence alignments, inflate assembly fragmentation, or introduce false variants.

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

Review the fastp JSON summary and the post-filtering FastQC reports to verify that low-quality tails have been successfully removed while ensuring an acceptable proportion of clean, usable read pairs remains for downstream applications.

---

## 3. Host Decontamination (`host_removal.sh`)

In clinical or host-associated metagenomics, a significant portion of raw reads often originate from the host organism rather than the target microbial or viral community; removing these host sequences reduces computational overhead and prevents misattribution during profiling. This step aligns clean reads against the human reference genome (`GRCh38_noalt_as`) using `Bowtie2` and extracts unmapped (host-free) paired reads via `samtools`.

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

Inspect the Bowtie2 alignment logs and MultiQC summary metrics to determine the exact percentage of host reads filtered out per sample, ensuring that a sufficient non-host fraction remains for assembly and profiling.

---

## 4. De Novo Assembly (`assembly.sh`)

De novo assembly reconstructs overlapping reads into longer contiguous sequences (contigs), allowing us to recover novel viral genomes, bacterial strains, or functional genes without relying strictly on existing reference databases.
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

Evaluate the assembly outputs by checking metrics like contig count, maximum length, and N50 to gauge the overall continuity and structural quality of the reconstructed genomic drafts.

---

## 5. Taxonomic Profiling (`kraken2.sh`)

Kraken2 allows for rapid, k-mer-based taxonomic classification of sequencing reads against a reference database, giving us an immediate quantitative overview of the microbial composition, relative abundances, and presence of specific pathogens within the sample.

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

Examine the Kraken2 report files to check the classification rates, verify the proportion of unclassified versus classified reads, and identify dominant taxa or specific pathogen signatures present in each sample.

---

## 6. Reference-Based Alignment (`mapping.sh`)
Reference-based alignment maps our host-free reads directly against a specific reference genome of interest, enabling deep variant calling, consensus generation, coverage depth assessment, and micro-diversity tracking for target organisms.

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

Inspect the generated sorted BAM files and alignment statistics to check genome-wide coverage breadth, mean sequencing depth, and mapping rates against the target reference accession.
