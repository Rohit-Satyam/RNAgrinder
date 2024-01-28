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
If you wanna run Kraken2 and Braken as well on your dataset to get a metagenomic blueprint, you can enable the Kraken2 classification as follows:

```bash
sbatch -N 1 -J sarasJob --mem=350G --time=3-24:00 --cpus-per-task=120 --mail-user=rohit.satyam@kaust.edu.sa --mail-type=FAIL \
--partition=batch -o sara.out -e sara.err --wrap="nextflow run main.nf --input 'data/*_L001_R{1,2}_001.fastq.gz' --outdir results \
 --mode PE --cpus 120 --k2db /ibex/scratch/projects/c2077/rohit/backup_runs/RNAgrinder/kraken2/index \
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
