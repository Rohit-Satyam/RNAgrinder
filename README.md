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
