# PROJECT: Processing Exome dataset of 3 samples with only chromsome 10

---

- Table of contents
  - Table of contents
  - Upstream analysis
    - Prepare
      - Setup working dir
      - Reference
      - Download tools
      - Download data
    - Raw data processing (FastQC & Trimmomatic)
    - Mapping (bwa)
    - Mapped post-processing & Alignment QC (GATK)
  - Downstream analysis
    - Variant calling (GATK - Haplotypecaller)
    - Annotation (Annovar)

## SETUP WORKING DIR

```bash
mkdir sra
mkdir raw
mkdir -p raw/qc_check
mkdir trimmed
mkdir -p trimmed/qc_check
mkdir bam
mkdir annotation
mkdir ref
mkdir -p ref/annotation
mkdir -p ref/genome
```

### Results

```text
.
├── align
├── raw
│   └── qc_check
├── ref
│   ├── annotation
│   └── genome
├── sra
└── trim
    └── qc_check
```
## DOWNLOAD REFFERENCE

### hg38 human genome fasta

```bash
wget  https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.fullAnalysisSet.chroms.tar.gz
tar xvzf hg38.fullAnalysisSet.chroms.tar.gz
```


## DOWNLOAD DATA

### Create a list of SRA Accession Number

```bash
nano SraAccList.csv
```

### Then copy and paste these Acessions to your SraAccList.csv

```text
SRR925745
SRR925746
SRR925747
```

## COMMAND TO DOWNLOAD RAW DATA OF THE PROJECT

```bash
for sra_acc in $(cat SraAccList.csv); do

 # Get the SRA
 prefetch -v $sra_acc --output-directory sra
    echo "=== $sra_acc sra file downloaded ==="
 # Download fastq file
 fastq-dump \
    --outdir raw/ \
    --split-files sra/${sra_acc}/${sra_acc}.sra \
    && gzip -f raw/*.fastq
    echo "=== Split $sra_acc file completed ==="

done
```

### Result

```bash
tree -h raw/

├── [4.0K]  qc_check
├── [4.77G]  SRR925745_1.fastq.gz
├── [4.57G]  SRR925745_2.fastq.gz
├── [4.16G]  SRR925746_1.fastq.gz
├── [4.11G]  SRR925746_2.fastq.gz
├── [4.49G]  SRR925747_1.fastq.gz
├── [4.32G]  SRR925747_2.fastq.gz
```

## RAW DATA PROCESSING: QUALITY CONTROL

Use FastQC

```BASH
for file in $(ls raw/*gz); do
    fastqc $file -o raw/qc_check/
done
```

## RAW DATA PROCESSING: TRIMMING & FILTERING

Trim and remove bad quality sequence with Trimmomatic

```bash
for file in $(ls raw/*_1.fastq.gz);do
    # Extract the base name without the extension
    prefix="raw/"
    foo=${file#"$prefix"}
    base_name=${foo%_*.*.*}

    # Construct the input and output file names
    read_1="${base_name}_1.fastq.gz"
    read_2="${base_name}_2.fastq.gz"
    output_paired_1="${base_name}_R1_paired.fastq.gz"
    output_unpaired_1="${base_name}_R1_unpaired.fastq.gz"
    output_paired_2="${base_name}_R2_paired.fastq.gz"
    output_unpaired_2="${base_name}_R2_unpaired.fastq.gz"
    java -jar $trimmomatic_jar PE \
      raw/$read_1 raw/$read_2 \
      trim/$output_paired_1 trim/$output_unpaired_1 \
       trim/$output_paired_2 trim/$output_unpaired_2 \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 \
    
    # Remove unpaired file
    rm trim/$output_unpaired_1 trim/$output_unpaired_2

    # Check fastqc again
    fastqc trim/$output_paired_1 -o $p_trim/qc_check
    fastqc trim/$output_paired_2 -o $p_trim/qc_check
done
```

## ALIGN

### Index the genome

```bash
#index only chromosome 10 fasta file
bwa index -a bwtsw ref/chr10.fa
```

### Alignment

```bash
bwa mem $p_ref/chr10.fa \
trim/SRR925745_R1_paired.fastq.gz \
trim/SRR925745_R2_paired.fastq.gz |\
samtools view -bS - \
samtools sort - \
samtools view -C -T ref/chr10.fa \
-o bam/aln-pe745.cram -
```


## POST-PROCESSING ALIGNMENT DATA: REMOVE DUPLICATE

```BASH
#sort file
samtools view -bu -T ref/chr10.fa \
bam/aln-pe745.cram |\
java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar SortSam \
  -I dev/stdin \
  --OUTPUT bam/sorted-pe745.bam \
  --SORT_ORDER coordinate
# keep only read mapped to chromosome 10
samtools view -bo bam/sortedchr10-pe745.bam bam/sorted-pe745.bam chr10 
#add Readgroup
picard-tools AddOrReplaceReadGroups \
       I=bam/sortedchr10-pe745.bam \
       O=bam/sortedchr10-RGpe745.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20
#mark duplicate
java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar MarkDuplicates \
  --INPUT bam/sortedchr10-RGpe745.bam.bam \
  --OUTPUT bam/markdup-sortedchr10-pe745.bam \
  --METRICS_FILE bam/markdup-sortedchr10-pe745.txt
#remove dup
java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar MarkDuplicates \
  --INPUT bam/sortedchr10-RGpe745.bam.bam \
  --OUTPUT bam/rmdup-sortedchr10-pe745.bam \
  --METRICS_FILE bam/rmdup-sortedchr10-pe745.txt
--REMOVE_DUPLICATES true
```

## POST-PROCESSING ALIGNMENT DATA: QUALITY CONTROL 

    ```BASH
    #create dict file for refference
    java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar CreateSequenceDictionary \
       -R ref/chr10.fa \
       -O ref/chr10.dict
    #1st baserecalibrator => recal.table 
    java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar BaseRecalibrator \
      -R ref/chr10.fa \
      -I bam/rmdup-sortedchr10-pe745.bam \
      -known-sites ref/Hg38.dbsnp138.vcf \
      -O bam/beforrecal-pe745.table 
    #apply bqsr
    java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar ApplyBQSR \
      -R ref/chr10.fa \
      -I bam/rmdup-sortedchr10-pe745.bam \
      -bqsr bam/beforrecal-pe745.table \
      -O bam/bqsr-rmdup-sortedchr10-pe745.bam
    #2nd baserecalibrator recal bam file
    java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar BaseRecalibrator \
      -R ref/chr10.fa \
      -I bam/bqsr-rmdup-sortedchr10-pe745.bam \
      -known-sites ref/Hg38.dbsnp138.vcf \
      -O bam/afterrecal-pe745.table
    java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar AnalyzeCovariates \
      -before bam/beforrecal-pe745.table  \
      -after bam/afterrecal-pe745.table  \
      -plots bam/recal-pe745-plots.pdf
      
    ```

### CALCULATE COVERAGE


## DOWNSTREAM ANALYSIS

### CALLING VARIANT

```bash
java -jar gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar HaplotypeCaller \
  -R ref/chr10.fa  \
  -I bam/bqsr-rmdup-sortedchr10-pe745.bam \
  -O pe745.chr10hg38.vcf.gz
```
### ANNOTATION

```bash
#download appropriate databases
annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
annovar/annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 
annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20211019 humandb/
annovar/table_annovar.pl pe745.chr10hg38.vcf \
humandb/ -buildver hg38 \
-out pe745.chr10hg38-clin.vcf \
-remove -protocol clinvar_20211019 \
-operation f -nastring . -vcfinput
annovar/table_annovar.pl pe745.chr10hg38-clin.vcf.avinput \
humandb/ -buildver hg38 \
-out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
-operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
```
