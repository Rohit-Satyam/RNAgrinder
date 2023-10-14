# RNAgrinder
A nextflow workflow for RNA-Seq Data assessment


# Preparing the Input Files 
Usually no preparation is required for running the pipeline but if you want to index your genome using STAR custom indexing parameters like `--sjdbOverhang`, you 
can do do and provide the index directory path using `--index_dir`.

**Note for non GENCODE GTF users** 
If using GTF file for any other organism other than Human, please check if the `transcript_type` and `gene_type` tags are present in the GTF. Though these tags are not required by STAR,
they are required by `gtexCollapseAnnotation.py` to alter GTF file used by RNASeQC. For example, in Plasmodium GTF file (GTF file produced from GFF using `AGAT`), the above mentioned 
tags are usualy missing and can be added as follows (as per discussion [here](https://github.com/NBISweden/AGAT/issues/398)

```
sed -i 's/gene_ebi_biotype/transcript_type/g'  PlasmoDB-64_Pfalciparum3D7.gtf
sed -i 's/ebi_biotype/gene_type/g'  PlasmoDB-64_Pfalciparum3D7.gtf
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
--gtf /home/subudhak/Documents/amit_timeseries_redo/reference/PlasmoDB-64_Pfalciparum3D7.gtf \
--mode SE --rrnaUse ribodetector --index_dir 00_index/
```

## How to run multiple samples on IBEX using pipeline
Submitting all samples using nextflow using regex `*R{1,2}*` will still be time consuming no matter how many worker threads you use. So I use divide and conquer strategy. I fire a separate job for each samples as followd

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
