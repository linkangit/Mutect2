# Complete Beginner's Guide to Somatic Variant Calling with GATK4 Mutect2

Welcome to the comprehensive guide for identifying cancer-causing mutations in tumor samples! This tutorial will walk you through the entire process of somatic variant calling - from raw sequencing data to annotated mutations ready for biological interpretation.

## 🧬 What Are We Doing?

**Somatic variant calling** is the process of identifying genetic mutations that occur specifically in cancer cells (tumor) but not in normal healthy cells. These mutations can drive cancer development and are crucial for:

- Understanding cancer biology
- Developing targeted therapies
- Personalizing treatment plans
- Research into cancer mechanisms

We'll be comparing DNA from **tumor tissue** against **normal tissue** from the same patient to find mutations that are unique to the cancer.

## 📚 Background Knowledge

### What is GATK4 Mutect2?
GATK (Genome Analysis Toolkit) is the gold standard for variant calling in genomics. Mutect2 is specifically designed for finding somatic mutations by comparing tumor and normal samples.

### Key Concepts You Should Know:
- **FASTQ files**: Raw sequencing data containing DNA sequences and quality scores
- **Reference genome**: The "standard" human genome we compare our samples against
- **BAM files**: Binary alignment files showing where each DNA read maps to the reference
- **VCF files**: Variant Call Format files containing information about genetic variants
- **Read alignment**: The process of matching sequencing reads to their location in the reference genome

## 🛠️ Prerequisites

### Required Software
Before starting, you'll need these tools installed:

```bash
# Core alignment and processing tools
sudo apt-get update
sudo apt-get install -y bwa samtools

# Quality control
sudo apt-get install -y fastqc

# Read processing
sudo apt-get install -y trimmomatic

# GATK4 (download from Broad Institute)
wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
unzip gatk-4.4.0.0.zip
# Add to PATH or use full path to gatk executable

# seqtk for read subsampling
sudo apt-get install -y seqtk
```

### System Requirements
- **RAM**: Minimum 16GB, recommended 32GB+
- **Storage**: At least 100GB free space for this tutorial
- **CPU**: Multi-core processor recommended (4+ cores)
- **OS**: Linux or macOS (Windows users can use WSL)

### Computational Time Expectations
- **Full workflow with subset data**: 2-4 hours
- **Full workflow with complete WGS**: 12-24 hours
- **Individual steps**: 10 minutes to 2 hours each

## 📁 Understanding Our Data

We're using real cancer genomics data from the Genome in a Bottle (GIAB) consortium:

- **Normal sample**: `HG008-N-D` (healthy duodenal tissue)
- **Tumor sample**: `HG008-T` (pancreatic cancer cell line)
- **Sequencing type**: Whole genome sequencing (WGS)
- **Read length**: 150bp paired-end
- **Coverage**: ~30x (each position covered ~30 times on average)

## 🚀 Step-by-Step Tutorial

### Step 1: Project Setup and Data Download

Let's start by creating our project structure and downloading the sequencing data.

```bash
# Create our main project directory
mkdir -p Project/reads
cd Project/reads

echo "📥 Downloading sequencing data - this will take 15-30 minutes..."

# Normal tissue sample (healthy control)
echo "Downloading normal sample R1 (forward reads)..."
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz

echo "Downloading normal sample R2 (reverse reads)..."
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz

# Tumor tissue sample (cancer)
echo "Downloading tumor sample R1..."
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz

echo "Downloading tumor sample R2..."
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/NYGC_Illumina-WGS_20231023/HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz

cd ..
echo "✅ Data download complete!"
```

**What's happening here?**
- We're downloading paired-end reads (R1 and R2) for both samples
- Each file contains millions of DNA sequences from the sequencing machine
- The normal sample will be our baseline for comparison
- The tumor sample is where we expect to find cancer mutations

### Step 2: Create Test Dataset (Recommended for Beginners)

Working with full genome data can be overwhelming and slow. Let's create smaller test files to learn the workflow:

```bash
cd Project/reads

echo "🔬 Creating smaller test datasets for faster processing..."

# Create subset with 1 million reads each (much faster for learning)
seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R1.fastq.gz 1000000 > HG008-N-D_subset_R1.fastq
seqtk sample -s100 HG008-N-D_CGGACAAC-AATCCGGA_H3LLJDSXC_L001_001.R2.fastq.gz 1000000 > HG008-N-D_subset_R2.fastq
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R1.fastq.gz 1000000 > HG008-T_subset_R1.fastq
seqtk sample -s100 HG008-T_TTCCTGTT-AAGATACT_HJVY2DSX7_L001_001.R2.fastq.gz 1000000 > HG008-T_subset_R2.fastq

cd ..
echo "✅ Test datasets created! Using these will make the tutorial much faster."
```

**Why subset the data?**
- Full WGS files are 50-100GB each and take hours to process
- 1M reads gives us enough data to learn without waiting all day
- You can always repeat with full data once you understand the process

### Step 3: Prepare the Reference Genome

The reference genome is like a "map" that we'll use to figure out where each sequencing read belongs.

```bash
mkdir -p Project/supporting_files/hg38
cd Project/supporting_files/hg38

echo "📚 Downloading human reference genome (hg38)..."
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

echo "🔍 Creating reference genome index files..."
# These index files help tools quickly find specific genome regions
samtools faidx hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

echo "📋 Downloading known variant sites for quality control..."
# These are previously discovered variations in human populations
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

cd ../..
echo "✅ Reference genome preparation complete!"
```

**What are these files?**
- **hg38.fa**: The human reference genome sequence
- **hg38.fa.fai**: Index file for fast random access
- **hg38.dict**: Dictionary of chromosome names and lengths
- **dbsnp138.vcf**: Database of known genetic variations

### Step 4: Quality Control Check

Let's examine the quality of our sequencing data before processing:

```bash
mkdir -p Project/reads/fastqc

echo "🔍 Running quality control analysis..."
fastqc Project/reads/HG008-N-D_subset_R1.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-N-D_subset_R2.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-T_subset_R1.fastq -o Project/reads/fastqc
fastqc Project/reads/HG008-T_subset_R2.fastq -o Project/reads/fastqc

echo "✅ Quality control complete! Check the HTML reports in Project/reads/fastqc/"
echo "📊 Look for: sequence quality scores, GC content, adapter contamination"
```

**Understanding FastQC Reports:**
- **Per base sequence quality**: Should be green (good) across most of the read length
- **Adapter content**: Should be minimal (blue line near bottom)
- **GC content**: Should follow expected distribution for human DNA

### Step 5: Clean Up the Reads

Real sequencing data contains artifacts that can interfere with analysis. Let's clean it up:

```bash
mkdir -p Project/reads/trimmed

echo "📦 Downloading adapter sequences..."
wget -O Project/reads/TruSeq3-PE.fa https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

echo "✂️ Trimming adapters and low-quality bases from normal sample..."
trimmomatic PE -threads 4 \
  Project/reads/HG008-N-D_subset_R1.fastq Project/reads/HG008-N-D_subset_R2.fastq \
  Project/reads/trimmed/HG008-N-D_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-N-D_trimmed_R1_unpaired.fastq \
  Project/reads/trimmed/HG008-N-D_trimmed_R2_paired.fastq Project/reads/trimmed/HG008-N-D_trimmed_R2_unpaired.fastq \
  ILLUMINACLIP:Project/reads/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "✂️ Trimming adapters and low-quality bases from tumor sample..."
trimmomatic PE -threads 4 \
  Project/reads/HG008-T_subset_R1.fastq Project/reads/HG008-T_subset_R2.fastq \
  Project/reads/trimmed/HG008-T_trimmed_R1_paired.fastq Project/reads/trimmed/HG008-T_trimmed_R1_unpaired.fastq \
  Project/reads/trimmed/HG008-T_trimmed_R2_paired.fastq Project/reads/trimmed/HG008-T_trimmed_R2_unpaired.fastq \
  ILLUMINACLIP:Project/reads/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "✅ Read trimming complete!"
```

**What does trimming do?**
- Removes adapter sequences (artificial DNA added during sequencing)
- Cuts off low-quality ends of reads
- Filters out reads that become too short after trimming
- Improves alignment accuracy and reduces false positive variants

### Step 6: Align Reads to Reference Genome

Now we'll map each sequencing read to its location in the human genome:

```bash
mkdir -p Project/aligned

echo "🗺️ Creating BWA index for the reference genome..."
bwa index Project/supporting_files/hg38/hg38.fa

echo "🎯 Aligning normal sample reads to reference genome..."
bwa mem -t 4 -R "@RG\tID:HG008-N-D\tPL:ILLUMINA\tSM:HG008-N-D" \
  Project/supporting_files/hg38/hg38.fa \
  Project/reads/trimmed/HG008-N-D_trimmed_R1_paired.fastq \
  Project/reads/trimmed/HG008-N-D_trimmed_R2_paired.fastq | \
  samtools view -bS - | samtools sort -o Project/aligned/HG008-N-D.sorted.bam

echo "🎯 Aligning tumor sample reads to reference genome..."
bwa mem -t 4 -R "@RG\tID:HG008-T\tPL:ILLUMINA\tSM:HG008-T" \
  Project/supporting_files/hg38/hg38.fa \
  Project/reads/trimmed/HG008-T_trimmed_R1_paired.fastq \
  Project/reads/trimmed/HG008-T_trimmed_R2_paired.fastq | \
  samtools view -bS - | samtools sort -o Project/aligned/HG008-T.sorted.bam

echo "📇 Creating index files for quick access..."
samtools index Project/aligned/HG008-N-D.sorted.bam
samtools index Project/aligned/HG008-T.sorted.bam

echo "✅ Read alignment complete!"
```

**What's happening in this step?**
- **BWA-MEM**: Aligns reads to find their most likely genomic location
- **Read groups**: Add metadata identifying each sample
- **Sorting**: Orders reads by genomic position for efficient processing
- **Indexing**: Creates index files for fast random access to any genome region

### Step 7: Mark Duplicate Reads

Some reads are sequenced multiple times (PCR duplicates). We need to mark these to avoid counting them multiple times:

```bash
echo "🔍 Marking duplicate reads in normal sample..."
gatk MarkDuplicatesSpark \
  -I Project/aligned/HG008-N-D.sorted.bam \
  -O Project/aligned/HG008-N-D_sorted_dedup_reads.bam

echo "🔍 Marking duplicate reads in tumor sample..."
gatk MarkDuplicatesSpark \
  -I Project/aligned/HG008-T.sorted.bam \
  -O Project/aligned/HG008-T_sorted_dedup_reads.bam

echo "📇 Indexing deduplicated files..."
samtools index Project/aligned/HG008-N-D_sorted_dedup_reads.bam
samtools index Project/aligned/HG008-T_sorted_dedup_reads.bam

echo "✅ Duplicate marking complete!"
```

**Why mark duplicates?**
- PCR amplification can create multiple copies of the same DNA fragment
- These would appear as multiple reads supporting the same variant
- Marking (not removing) prevents over-counting while preserving coverage information

### Step 8: Base Quality Score Recalibration

Sequencing machines sometimes systematically over- or under-estimate quality scores. Let's fix this:

```bash
echo "📊 Building quality recalibration model for normal sample..."
gatk BaseRecalibrator \
  -I Project/aligned/HG008-N-D_sorted_dedup_reads.bam \
  -R Project/supporting_files/hg38/hg38.fa \
  --known-sites Project/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  -O Project/aligned/HG008-N-D_recal_data.table

echo "📊 Building quality recalibration model for tumor sample..."
gatk BaseRecalibrator \
  -I Project/aligned/HG008-T_sorted_dedup_reads.bam \
  -R Project/supporting_files/hg38/hg38.fa \
  --known-sites Project/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  -O Project/aligned/HG008-T_recal_data.table

echo "🔧 Applying quality recalibration to normal sample..."
gatk ApplyBQSR \
  -I Project/aligned/HG008-N-D_sorted_dedup_reads.bam \
  -R Project/supporting_files/hg38/hg38.fa \
  --bqsr-recal-file Project/aligned/HG008-N-D_recal_data.table \
  -O Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam

echo "🔧 Applying quality recalibration to tumor sample..."
gatk ApplyBQSR \
  -I Project/aligned/HG008-T_sorted_dedup_reads.bam \
  -R Project/supporting_files/hg38/hg38.fa \
  --bqsr-recal-file Project/aligned/HG008-T_recal_data.table \
  -O Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam

echo "📇 Final indexing..."
samtools index Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam
samtools index Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam

echo "✅ Base quality recalibration complete!"
```

**What does BQSR accomplish?**
- Corrects systematic errors in quality score reporting
- Uses known variant sites to model true vs. false variants
- Improves accuracy of downstream variant calling
- Essential for high-quality variant detection

### Step 9: Quality Assessment

Let's check how well our alignment worked:

```bash
echo "📈 Collecting alignment statistics..."

gatk CollectAlignmentSummaryMetrics \
  R=Project/supporting_files/hg38/hg38.fa \
  I=Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  O=Project/aligned/HG008-N-D_alignment_metrics.txt

gatk CollectInsertSizeMetrics \
  INPUT=Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  OUTPUT=Project/aligned/HG008-N-D_insert_size_metrics.txt \
  HISTOGRAM_FILE=Project/aligned/HG008-N-D_insert_size_histogram.pdf

gatk CollectAlignmentSummaryMetrics \
  R=Project/supporting_files/hg38/hg38.fa \
  I=Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  O=Project/aligned/HG008-T_alignment_metrics.txt

gatk CollectInsertSizeMetrics \
  INPUT=Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  OUTPUT=Project/aligned/HG008-T_insert_size_metrics.txt \
  HISTOGRAM_FILE=Project/aligned/HG008-T_insert_size_histogram.pdf

echo "✅ Quality assessment complete! Check the metrics files for alignment statistics."
```

### Step 10: Download Resources for Variant Calling

Mutect2 needs additional resources to distinguish real somatic mutations from artifacts:

```bash
mkdir -p Project/supporting_files/mutect2_supporting_files
cd Project/supporting_files/mutect2_supporting_files

echo "📚 Downloading gnomAD population frequency data..."
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

echo "🛡️ Downloading Panel of Normals..."
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

echo "🎯 Downloading exome intervals..."
wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list

cd ../..
echo "✅ Variant calling resources downloaded!"
```

**What are these resources?**
- **gnomAD**: Population allele frequencies to filter common germline variants
- **Panel of Normals**: Common artifacts seen across many normal samples
- **Exome intervals**: Regions of the genome where we expect to find variants

### Step 11: Set Up Variant Annotation

We'll need Funcotator to annotate our variants with biological information:

```bash
mkdir -p Project/supporting_files/funcotator
cd Project/supporting_files/funcotator

echo "📚 Downloading Funcotator annotation databases..."
wget https://console.cloud.google.com/storage/browser/_details/broad-public-datasets/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
tar -xzf funcotator_dataSources.v1.8.hg38.20230908s.tar.gz

cd ../..
echo "✅ Annotation resources ready!"
```

### Step 12: The Main Event - Call Somatic Variants!

This is where we actually identify mutations that are present in tumor but not normal tissue:

```bash
mkdir -p Project/results

echo "🔬 Running Mutect2 somatic variant calling..."
echo "This step identifies mutations unique to the tumor sample..."

gatk Mutect2 \
  -R Project/supporting_files/hg38/hg38.fa \
  -I Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  -I Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  -tumor HG008-T \
  -normal HG008-N-D \
  --germline-resource Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals Project/supporting_files/mutect2_supporting_files/1000g_pon.hg38.vcf.gz \
  -O Project/results/HG008_somatic_variants_mutect2.vcf.gz \
  --f1r2-tar-gz Project/results/HG008_f1r2.tar.gz

echo "✅ Initial variant calling complete!"
```

**What is Mutect2 doing?**
- Compares tumor and normal samples at every genomic position
- Uses statistical models to determine if differences are real mutations
- Considers read depth, base quality, and population frequencies
- Outputs candidate somatic variants for further filtering

### Step 13: Estimate Sample Cross-Contamination

Sometimes tumor and normal samples can be contaminated with each other. Let's check:

```bash
echo "🔍 Checking for cross-sample contamination..."

echo "Analyzing tumor sample..."
gatk GetPileupSummaries \
  -I Project/aligned/HG008-T_sorted_dedup_bqsr_reads.bam \
  -V Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  -L Project/supporting_files/mutect2_supporting_files/exome_calling_regions.v1.1.interval_list \
  -O Project/results/HG008-T_getpileupsummaries.table

echo "Analyzing normal sample..."
gatk GetPileupSummaries \
  -I Project/aligned/HG008-N-D_sorted_dedup_bqsr_reads.bam \
  -V Project/supporting_files/mutect2_supporting_files/af-only-gnomad.hg38.vcf.gz \
  -L Project/supporting_files/mutect2_supporting_files/exome_calling_regions.v1.1.interval_list \
  -O Project/results/HG008-N-D_getpileupsummaries.table

echo "Calculating contamination levels..."
gatk CalculateContamination \
  -I Project/results/HG008-T_getpileupsummaries.table \
  -matched Project/results/HG008-N-D_getpileupsummaries.table \
  -O Project/results/HG008_pair_calculatecontamination.table

echo "✅ Contamination analysis complete!"
```

### Step 14: Model Sequencing Artifacts

Some mutations might actually be artifacts from the sequencing process. Let's model these:

```bash
echo "🔧 Learning read orientation artifacts..."
gatk LearnReadOrientationModel \
  -I Project/results/HG008_f1r2.tar.gz \
  -O Project/results/read-orientation-model.tar.gz

echo "✅ Artifact modeling complete!"
```

### Step 15: Filter the Variant Calls

Now we'll apply all our quality filters to get high-confidence somatic variants:

```bash
echo "🚰 Filtering variants to remove false positives..."
gatk FilterMutectCalls \
  -V Project/results/HG008_somatic_variants_mutect2.vcf.gz \
  -R Project/supporting_files/hg38/hg38.fa \
  --contamination-table Project/results/HG008_pair_calculatecontamination.table \
  --ob-priors Project/results/read-orientation-model.tar.gz \
  -O Project/results/HG008_somatic_variants_filtered_mutect2.vcf

echo "✅ Variant filtering complete!"
```

**What gets filtered out?**
- Variants with low quality scores
- Variants that appear to be sequencing artifacts
- Variants likely due to sample contamination
- Variants that don't meet statistical significance thresholds

### Step 16: Annotate the Final Variants

Finally, let's add biological information to our variants:

```bash
echo "📝 Adding biological annotations to variants..."
gatk Funcotator \
  --variant Project/results/HG008_somatic_variants_filtered_mutect2.vcf \
  --reference Project/supporting_files/hg38/hg38.fa \
  --ref-version hg38 \
  --data-sources-path Project/supporting_files/funcotator/funcotator_dataSources.v1.8.hg38.20230908s \
  --output Project/results/HG008_somatic_variants_functotated.vcf \
  --output-file-format VCF

echo "🎉 WORKFLOW COMPLETE! 🎉"
echo ""
echo "Your results are in Project/results/:"
echo "  - HG008_somatic_variants_filtered_mutect2.vcf (filtered variants)"
echo "  - HG008_somatic_variants_functotated.vcf (annotated variants)"
```

## 📊 Understanding Your Results

### Key Output Files

1. **HG008_somatic_variants_filtered_mutect2.vcf**
   - High-confidence somatic variants
   - Includes quality metrics and filter status
   - Ready for downstream analysis

2. **HG008_somatic_variants_functotated.vcf**
   - Same variants with biological annotations
   - Gene names, protein effects, clinical significance
   - Ready for interpretation

### Interpreting VCF Files

VCF files contain these key columns:
- **CHROM**: Chromosome where variant is located
- **POS**: Position on the chromosome
- **REF**: Reference genome sequence
- **ALT**: Alternative (mutated) sequence
- **QUAL**: Quality score for the variant
- **FILTER**: Whether the variant passed quality filters
- **INFO**: Additional information about the variant

### Example Results Analysis

```bash
# Count how many variants we found
echo "Total variants found:"
grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | wc -l

# Count variants that passed filters
echo "High-quality variants:"
grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | grep "PASS" | wc -l

# Look at the first few variants
echo "First 5 variants:"
grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | head -5
```

## 🔧 Troubleshooting Common Issues

### Memory Issues
```bash
# If you get out-of-memory errors, try reducing thread count:
# Change -threads 4 to -threads 2
# Or increase Java heap space:
export GATK_LOCAL_JAR="/path/to/gatk.jar"
export _JAVA_OPTIONS="-Xmx8g"  # Use 8GB RAM
```

### Slow Processing
```bash
# Process specific chromosomes only (much faster):
gatk Mutect2 \
  -L chr20 \  # Only analyze chromosome 20
  [other parameters...]
```

### File Permission Errors
```bash
# Make sure you have write permissions:
chmod -R 755 Project/
```

## 🎯 Next Steps

### Biological Interpretation
1. **Filter by impact**: Focus on variants in protein-coding regions
2. **Check known databases**: Search COSMIC, ClinVar for known cancer variants
3. **Pathway analysis**: See which biological pathways are affected
4. **Literature review**: Research specific genes with mutations

### Advanced Analyses
1. **Mutational signatures**: Identify patterns in your mutations
2. **Copy number analysis**: Look for chromosomal gains/losses
3. **Structural variants**: Find large rearrangements
4. **Tumor purity estimation**: Determine cancer cell fraction

### Validation
1. **PCR validation**: Confirm important variants with targeted sequencing
2. **Orthogonal methods**: Use different variant callers for comparison
3. **Functional studies**: Test variant effects in cell culture

## 📚 Additional Resources

### Learning More
- [GATK Documentation](https://gatk.broadinstitute.org/hc/en-us)
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
- [SAMtools Documentation](http://samtools.sourceforge.net/)
- [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

### Cancer Genomics Databases
- [COSMIC](https://cancer.sanger.ac.uk/cosmic) - Catalogue of Somatic Mutations in Cancer
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) - Clinical significance of variants
- [gnomAD](https://gnomad.broadinstitute.org/) - Population allele frequencies
- [cBioPortal](https://www.cbioportal.org/) - Cancer genomics data portal
- [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) - The Cancer Genome Atlas

### Visualization Tools
- [IGV](https://software.broadinstitute.org/software/igv/) - Integrative Genomics Viewer
- [UCSC Genome Browser](https://genome.ucsc.edu/) - Web-based genome browser
- [Ensembl](https://www.ensembl.org/) - Genome annotation and analysis

## 💡 Pro Tips for Success

### Data Management
```bash
# Always organize your data clearly
Project/
├── reads/           # Raw sequencing data
├── supporting_files/ # Reference genomes and databases
├── aligned/         # Processed BAM files
├── results/         # Final variant calls
└── logs/           # Keep logs of all commands
```

### Keep Detailed Logs
```bash
# Log all your commands with timestamps
echo "$(date): Starting BWA alignment" >> Project/logs/analysis.log
bwa mem [parameters...] 2>&1 | tee -a Project/logs/bwa.log
```

### Version Control
```bash
# Document software versions
echo "Software versions used:" > Project/software_versions.txt
gatk --version >> Project/software_versions.txt
bwa 2>&1 | head -3 >> Project/software_versions.txt
samtools --version >> Project/software_versions.txt
```

### Backup Critical Files
```bash
# Always backup your results
cp -r Project/results/ Project/results_backup_$(date +%Y%m%d)/
```

## 🧪 Validation and Quality Control

### Manual Inspection of Variants
```bash
# Look at specific interesting variants in IGV
# Load your BAM files and VCF file to visually inspect variants
# Check for:
# - Adequate read depth (>10 reads)
# - Balanced strand support
# - No clustering of variants (potential artifacts)
# - Clean breakpoints for indels
```

### Statistical Validation
```bash
# Calculate basic statistics
echo "=== VARIANT CALLING SUMMARY ===" > Project/results/summary_stats.txt
echo "Date: $(date)" >> Project/results/summary_stats.txt
echo "" >> Project/results/summary_stats.txt

# Total variants called
total_vars=$(grep -v "^#" Project/results/HG008_somatic_variants_mutect2.vcf.gz | wc -l)
echo "Total variants called: $total_vars" >> Project/results/summary_stats.txt

# Variants that passed filtering
passed_vars=$(zgrep -v "^#" Project/results/HG008_somatic_variants_mutect2.vcf.gz | grep "PASS" | wc -l)
echo "High-quality variants (PASS): $passed_vars" >> Project/results/summary_stats.txt

# Calculate filtering rate
filter_rate=$(echo "scale=2; ($total_vars - $passed_vars) * 100 / $total_vars" | bc)
echo "Filtering rate: ${filter_rate}%" >> Project/results/summary_stats.txt
```

## 🔍 Advanced Analysis Examples

### Extract High-Impact Variants
```bash
# Find variants in protein-coding regions
grep -v "^#" Project/results/HG008_somatic_variants_functotated.vcf | \
grep "PASS" | \
grep -E "(MISSENSE|NONSENSE|FRAME_SHIFT|SPLICE_SITE)" > \
Project/results/high_impact_variants.txt

echo "High-impact variants saved to high_impact_variants.txt"
```

### Create Variant Summary Table
```bash
# Create a simple summary of variants by type
echo "Creating variant type summary..."
echo "Variant_Type,Count" > Project/results/variant_types.csv

# Count SNPs (single nucleotide changes)
snp_count=$(grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | \
           grep "PASS" | \
           awk 'length($4)==1 && length($5)==1' | wc -l)
echo "SNP,$snp_count" >> Project/results/variant_types.csv

# Count insertions
ins_count=$(grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | \
           grep "PASS" | \
           awk 'length($4)<length($5)' | wc -l)
echo "Insertion,$ins_count" >> Project/results/variant_types.csv

# Count deletions
del_count=$(grep -v "^#" Project/results/HG008_somatic_variants_filtered_mutect2.vcf | \
           grep "PASS" | \
           awk 'length($4)>length($5)' | wc -l)
echo "Deletion,$del_count" >> Project/results/variant_types.csv

echo "✅ Variant summary saved to variant_types.csv"
```

## ⚠️ Important Considerations

### Ethical and Legal Aspects
- **Patient consent**: Ensure proper consent for genomic analysis
- **Data privacy**: Follow institutional guidelines for handling genetic data
- **HIPAA compliance**: Protect patient information appropriately
- **Data sharing**: Understand restrictions on sharing genomic data

### Clinical Interpretation Warnings
- **This is for research/educational purposes only**
- **Clinical decisions should never be based solely on these results**
- **Always consult with qualified geneticists/oncologists**
- **Validate important findings with clinical-grade testing**

### Technical Limitations
- **Coverage bias**: Some genome regions are difficult to sequence
- **Repetitive regions**: May have alignment artifacts
- **Copy number variants**: Not detected by this workflow
- **Structural variants**: Large rearrangements not captured
- **Tumor heterogeneity**: Mixed populations of cancer cells

## 📈 Performance Expectations

### Resource Usage (1M read subset)
- **Runtime**: 2-4 hours on modern hardware
- **Peak RAM**: ~8-16GB
- **Disk space**: ~20-30GB temporary files
- **CPU cores**: Scales well with 4-8 cores

### Resource Usage (Full WGS)
- **Runtime**: 12-24 hours
- **Peak RAM**: ~32-64GB for large files
- **Disk space**: ~500GB-1TB temporary files
- **CPU cores**: Benefits from 16+ cores

### Expected Results (typical cancer sample)
- **Total variants**: 10,000-100,000 (before filtering)
- **High-quality variants**: 100-10,000 (after filtering)
- **Protein-affecting variants**: 50-500
- **High-impact variants**: 5-50

## 🤝 Getting Help

### Common Questions
1. **"My variant calling is taking forever"**
   - Use subset data first
   - Try single chromosome analysis (-L chr20)
   - Check system resources

2. **"I'm getting permission denied errors"**
   - Check file permissions with `ls -la`
   - Ensure write access to output directories

3. **"GATK is running out of memory"**
   - Increase Java heap: `export _JAVA_OPTIONS="-Xmx16g"`
   - Use fewer threads to reduce memory usage

4. **"I found too many/too few variants"**
   - Check input data quality
   - Verify normal/tumor sample labels
   - Review filtering parameters

### Where to Get Support
- **GATK Forum**: [gatk.broadinstitute.org/hc/en-us/community/topics](https://gatk.broadinstitute.org/hc/en-us/community/topics)
- **Bioinformatics Stack Exchange**: [bioinformatics.stackexchange.com](https://bioinformatics.stackexchange.com)
- **Local bioinformatics core**: Many institutions have bioinformatics support
- **GitHub Issues**: For software-specific problems

## 🎓 Learning Path Recommendations

### For Beginners
1. **Complete this tutorial** with subset data
2. **Learn VCF format** - understand your output files
3. **Try IGV** - visualize your variants
4. **Read GATK documentation** - understand the theory
5. **Repeat with full data** - apply what you've learned

### For Intermediate Users
1. **Try different variant callers** (MuSE, Strelka2, VarScan2)
2. **Learn variant annotation** (VEP, ANNOVAR)
3. **Explore mutational signatures** (SigProfiler, deconstructSigs)
4. **Add copy number analysis** (GATK CNV, PURPLE)
5. **Study cancer genomics literature**

### For Advanced Users
1. **Optimize pipeline performance** (parallelization, cloud computing)
2. **Implement quality control metrics** (contamination, tumor purity)
3. **Add structural variant calling** (Manta, Delly, GRIDSS)
4. **Integrate multi-omics data** (RNA-seq, ChIP-seq)
5. **Develop custom analysis scripts** (R, Python)

## 🌟 Final Thoughts

Congratulations! You've completed a comprehensive somatic variant calling workflow. This represents a significant achievement in computational biology. You now have the foundation to:

- Identify cancer-causing mutations from sequencing data
- Understand the quality control steps essential for accurate results  
- Interpret variant calling outputs and their biological significance
- Troubleshoot common issues in genomics pipelines

Remember that this is just the beginning. Cancer genomics is a rapidly evolving field with new tools, databases, and insights emerging regularly. Stay curious, keep learning, and don't hesitate to ask questions in the bioinformatics community.

**Happy variant hunting!** 🔬✨

---

## 📋 Quick Reference Commands

### Essential File Checks
```bash
# Check file sizes
ls -lh Project/*/

# Count reads in FASTQ
zcat file.fastq.gz | wc -l | awk '{print $1/4}'

# Check BAM file integrity
samtools quickcheck Project/aligned/*.bam

# Count variants in VCF
grep -v "^#" file.vcf | wc -l
```

### Pipeline Status Monitoring
```bash
# Monitor resource usage
htop
df -h  # Disk space
free -h  # Memory usage

# Check running processes
ps aux | grep gatk
ps aux | grep bwa
```

### Quick Quality Checks
```bash
# Check alignment rates
samtools flagstat Project/aligned/sample.bam

# View variant quality distribution
grep -v "^#" Project/results/*.vcf | cut -f6 | sort -n | tail -10
```

This tutorial provides a complete foundation for cancer genomics analysis. Practice with different datasets, explore the literature, and continue building your expertise in this exciting field!
