[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10518597.svg)](https://doi.org/10.5281/zenodo.10518597) <br>
If you use this pipeline, please cite it using

```
Rohit Satyam. (2024). Rohit-Satyam/RNAgrinder: v1.0.1 (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.10518597
```

![](pipeline.png)

A nextflow workflow for RNA-Seq Data assessment

# Procure the data
If you have your own data good. But if you don't and wanna fetch data from SRA, here is a way to do it quickly and efficiently using `fasterq-dump`, `parallel` and `pysradb`

```bash
## eg. First convert your GSE ID to SRP and then use it to get SRR ID of the entire project
pysradb gse-to-srp  GSE243125
pysradb srp-to-srr SRP460304 > meta.tsv

## now choose the SRR IDs you wanna download from meta.tsv and save them as srr.txt
## The use fasterq-dump and parallel
cat srr.txt | parallel -j 5 "fasterq-dump -e 20 --skip-technical --split-3 -p {}"
```

# Preparing the Input Files 
Usually, no preparation is required for running the pipeline but if you want to index your genome using STAR custom indexing parameters like `--sjdbOverhang`, you 
can do so and provide the index directory path using `--index_dir`. For GENCODE genomes, no editing of GTF files are required, but for non-model organisms, perform 
some GTF editing such as the example shown below:

**Note for non-GENCODE GTF users** 
If using the GTF file for any other organism other than Humans, please check if the `transcript_type` and `gene_type` tags are present in the GTF. Though these tags are not required by STAR,
they are required by `gtexCollapseAnnotation.py` to alter the GTF file used by RNASeQC. For example, in the Plasmodium GTF file (GTF file produced from GFF using `AGAT`), the above mentioned 
tags are usually missing and can be added as follows (as per discussion [here](https://github.com/NBISweden/AGAT/issues/398)

```bash
agat_convert_sp_gff2gtf.pl -gff PlasmoDB-66_Pfalciparum3D7.gff -o PlasmoDB-66_Pfalciparum3D7.gtf --gtf_version 3
sed -i 's/gene_ebi_biotype/transcript_type/g'  PlasmoDB-66_Pfalciparum3D7.gtf
sed -i 's/ebi_biotype/gene_type/g'  PlasmoDB-66_Pfalciparum3D7.gtf
```

# Adding ncRNA to your organism of interest

We will document here steps to convert the lncRNA annotation available publically for Plasmodium falciparum 3D7 and Toxoplasma gondii ME49 and how can be liftover these annotations using `Liftoff` tool.

### For Toxoplasma
```bash
# Because the annotation were generated from StringTie, the gene records were missing from GFF files. To add that, we will first use AGAT to fix this
agat_convert_sp_gxf2gxf.pl -g Toxolncrna_v59.gff -o fixed.gff
liftoff -g fixed.gff -o Toxolncrna_v68.gff  ToxoDB-68_TgondiiME49_Genome.fasta ToxoDB-59_TgondiiME49_Genome.fasta
## Replace "gene" and "transcript" with ncRNA_gene and lnc_RNA to match ToxoDB standards. Also, adding other descriptions
awk 'BEGIN{OFS="\t"}
$3 == "gene" { $9 = $9 ";description=lncRNA;ebi_biotype=lncRNA" } 
$3 == "transcript" { $9 = $9 ";description=lncRNA;gene_ebi_biotype=lncRNA" } 
{ print }' Toxolncrna_v68.gff > temp.gff

awk 'BEGIN{OFS="\t"} $3 == "gene" {$3 = "ncRNA_gene"} $3 == "transcript" {$3 = "lnc_RNA"} {print}' temp.gff

## Now let's fix the exon naming record since AGAT assigns a random name in the ID field if the exon name is missing, eg: agat-exon-1577. 
awk 'BEGIN{OFS="\t"} 
$3 == "exon" { 
    split($9, fields, ";"); 
    for (i in fields) {
        if (fields[i] ~ /^ID=/) {
            id_index = i;
        } else if (fields[i] ~ /^Parent=/) {
            split(fields[i], parent, "=");
            parent_value = parent[2];
        } else if (fields[i] ~ /^exon_number=/) {
            split(fields[i], exon_number, "=");
            exon_value = exon_number[2];
        }
    }
    if (id_index && parent_value && exon_value) {
        fields[id_index] = "ID=exon_" parent_value "-E" exon_value;
    }
    $9 = "";
    for (i = 1; i <= length(fields); i++) {
        if (fields[i] != "") {
            $9 = $9 ? $9 ";" fields[i] : fields[i];
        }
    }
}{ print }' revised.gff > ready2merge.gff

## Finally merging the two annotations and make a gtf file as well
agat_sp_merge_annotations.pl --gff ToxoDB-68_TgondiiME49.gff --gff ready2merge.gff --out ToxoDB-68_TgondiiME49_withlncRNA.gff
agat_convert_sp_gff2gtf.pl -gff ToxoDB-68_TgondiiME49_withlncRNA.gff -o ToxoDB-68_TgondiiME49_withlncRNA.gtf --gtf_version 3
```
Reference:
```
grep TGME49_500145 ToxoDB-68_TgondiiME49.gff
TGME49_chrVIIb	VEuPathDB	ncRNA_gene	1614480	1617748	.	+	.	ID=TGME49_500145;description=lncRNA;ebi_biotype=lncRNA
TGME49_chrVIIb	VEuPathDB	lnc_RNA	1614480	1617748	.	+	.	ID=TGME49_500145.R149;Parent=TGME49_500145;description=lncRNA;gene_ebi_biotype=lncRNA
TGME49_chrVIIb	VEuPathDB	exon	1614480	1617748	.	+	.	ID=exon_TGME49_500145.R149-E1;Parent=TGME49_500145.R149;gene_id=TGME49_500145
```
The ready2merge.gff looks like

```
TGME49_chrII	Liftoff	ncRNA_gene	7	2605	.	+	.	ID=MSTRG.1118;gene_id=MSTRG.1118;coverage=1.0;sequence_ID=1.0;extra_copy_number=0;copy_num_ID=MSTRG.1118_0;description=lncRNA;ebi_biotype=lncRNA
TGME49_chrII	Liftoff	lnc_RNA	7	2605	.	+	.	ID=MSTRG.1118.1;Parent=MSTRG.1118;gene_id=MSTRG.1118;transcript_id=MSTRG.1118.1;extra_copy_number=0;description=lncRNA;gene_ebi_biotype=lncRNA
TGME49_chrII	Liftoff	exon	7	623	.	+	.	ID=exon_MSTRG.1118.1-E1;Parent=MSTRG.1118.1;exon_number=1;gene_id=MSTRG.1118;transcript_id=MSTRG.1118.1;extra_copy_number=0
TGME49_chrII	Liftoff	exon	2298	2605	.	+	.	ID=exon_MSTRG.1118.1-E2;Parent=MSTRG.1118.1;exon_number=2;gene_id=MSTRG.1118;transcript_id=MSTRG.1118.1;extra_copy_number=0
TGME49_chrII	Liftoff	ncRNA_gene	825	11652	.	-	.	ID=MSTRG.1119;gene_id=MSTRG.1119;coverage=1.0;sequence_ID=1.0;extra_copy_number=0;copy_num_ID=MSTRG.1119_0;description=lncRNA;ebi_biotype=lncRNA
TGME49_chrII	Liftoff	lnc_RNA	825	11652	.	-	.	ID=MSTRG.1119.3;Parent=MSTRG.1119;gene_id=MSTRG.1119;transcript_id=MSTRG.1119.3;extra_copy_number=0;description=lncRNA;gene_ebi_biotype=lncRNA
TGME49_chrII	Liftoff	exon	825	1128	.	-	.	ID=exon_MSTRG.1119.3-E1;Parent=MSTRG.1119.3;exon_number=1;gene_id=MSTRG.1119;transcript_id=MSTRG.1119.3;extra_copy_number=0

```


## Command to run
For **PE Data**
```
nextflow run main.nf --input 'data/data_pe/*R{1,2}_001.fastq.gz' \
--outdir shuaibs_results --ref /home/subudhak/Documents/sara_samples/Iseq_COVID_batch2_and3_corrected_index/RNAgrinder/resources/GRCh38.primary_assembly.genome.fa \
--gtf /home/subudhak/Documents/sara_samples/Iseq_COVID_batch2_and3_corrected_index/RNAgrinder/resources/gencode.v43.primary_assembly.basic.annotation.gtf \
--mode PE --rrnaUse ribodetector --index_dir 00_index/
```
For **SE Data**

```
nextflow run main.nf --input 'data_se/*.fastq.gz'
--outdir results --ref /home/subudhak/Documents/amit_timeseries_redo/reference/PlasmoDB-64_Pfalciparum3D7_Genome.fasta \
--gtf /home/subudhak/Documents/amit_timeseries_redo/reference/PlasmoDB-66_Pfalciparum3D7.gtf \
--mode SE --rrnaUse ribodetector --index_dir 00_index/
```
If you want to run Kraken2 and Braken as well on your dataset to get a metagenomic blueprint, you can enable the Kraken2 classification as follows:

```bash
sbatch -N 1 -J sarasJob --mem=350G --time=3-24:00 --cpus-per-task=120 --mail-user=rohit.satyam@kaust.edu.sa --mail-type=FAIL \
--partition=batch -o sara.out -e sara.err --wrap="nextflow run main.nf --input 'data/*_L001_R{1,2}_001.fastq.gz' --outdir results \
 --mode PE --cpus 120 --skipKraken=false --k2db /ibex/scratch/projects/c2077/rohit/backup_runs/RNAgrinder/kraken2/index \
--ref $PWD/resources/PlasmoDB-66_Pfalciparum3D7_Genome.fasta --gtf $PWD/resources/PlasmoDB-66_Pfalciparum3D7.gtf"

```
However, you must first download the indexes of your interest (here we use PlusPF since that's the most comprehensive database) from their [webpage](https://benlangmead.github.io/aws-indexes/k2).
## How to run multiple samples on IBEX using a pipeline
Submitting all samples using nextflow using regex `*R{1,2}*` will still be time-consuming no matter how many worker threads you use. So I use the divide and conquer strategy. I fire a separate job for each sample as followed

**Step1**
```
## Create the input file list
ls -1 /ibex/scratch/projects/c2077/rohit/sara_novaseq/210801_A01018_0104_AHG7JCDSXY/Lane1/version_01/*R1*.gz > file
sed -i 's/_R1_/_*R{1,2}_/g' file
```

**Step2**
To fire multiple jobs we will use while loop

```
while read p
do
n=$(echo $p | xargs -n 1 basename | awk -F'_L001' '{print $1}')
sbatch -N 1 -J ${n}_lane1 --mem=100G --time=1-24:00 --cpus-per-task 12 --mail-user=rohit.satyam@kaust.edu.sa --mail-type=FAIL --partition=batch -o ${n}.out -e ${n}.err --wrap="nextflow run main.nf --input '$p' --outdir results_lane1 --ref /ibex/scratch/projects/c2077/rohit/RNAgrinder/resources/GRCh38.primary_assembly.genome.fa --gtf /ibex/scratch/projects/c2077/rohit/RNAgrinder/resources/gencode.v43.primary_assembly.basic.annotation.gtf --mode PE --rrnaUse ribodetector --index_dir /ibex/scratch/projects/c2077/rohit/RNAgrinder/resources/00_index/ --cpus 12 -w ${n}_work"
done < file
```

# To do list
1. Add library strandedness determination capability. Read https://groups.google.com/g/rna-star/c/mkooNLzyJYc and https://github.com/alexdobin/STAR/issues/1589 and try this tool: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04572-7
Tentative code
```
#!/bin/bash

# Loop through all ReadsPerGene.out.tab files in the current directory
for file in *.ReadsPerGene.out.tab; do
    echo "Processing $file..."
    # Calculate sums for each stranded column (assuming columns 2, 3, and 4 are unstranded, stranded-forward, and stranded-reverse, respectively)
    unstranded_sum=$(awk '{if(NR>4) sum+=$2} END {print sum}' $file)
    forward_stranded_sum=$(awk '{if(NR>4) sum+=$3} END {print sum}' $file)
    reverse_stranded_sum=$(awk '{if(NR>4) sum+=$4} END {print sum}' $file)
    
    # Compare the sums to infer strandedness
    if (( $(echo "$forward_stranded_sum > $reverse_stranded_sum" | bc -l) )); then
        strandedness="forward"
        column=3
    elif (( $(echo "$reverse_stranded_sum > $forward_stranded_sum" | bc -l) )); then
        strandedness="reverse"
        column=4
    else
        strandedness="unstranded"
        column=2
    fi
    
    # Output the decision
    echo "$file: Use column $column ($strandedness stranded)"
done

```
2. Add Gender recognition capability using https://rpubs.com/seungyeul/471026
