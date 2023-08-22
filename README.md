# RNAgrinder
A nextflow workflow for RNA-Seq Data assessment


# Preparing the Input Files 
### RSeQC rRNA BED file
RSeQC utility required rRNA BED file for Ribosomal RNA contamination. This can be obtained from [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) by setting the following values:
1. group: `All Tables`, database `hg38`, table `rmsk`.
2. Click on `filter` and add `rRNA` in front of `repClass` field.
3. For output choose either `BED-browser extensible data` (preferably) or `GTF`. I save it as `rrna.bed`
The BED file generated this way contains 6 columns. Column `chr`, `start` and `end` followed by `gene name`, `score` and `strand`. However, this BED file is not yet compatible. Add `start` and `end` columns one more time after `strand` filed. Next add two columns containing `0` and `1`. Then add one more column containing lengths of exons obtained by `end-start`. Now fill the last column with `0`. This gives us a standard 12 column BED format. For details about the columns refer to the document [here](https://agat.readthedocs.io/en/latest/gff_to_bed.html).

To achieve the above mentioned task using awk
 ```
 awk -F'\t' '{a=$3-$2; print $1,$2,$3,$4,$5,$6,$2,$3,0,1,a,0}' rrna.bed > rrna_mod.bed
 ```
> NOTE: I tried using AGAT `agat_convert_sp_gff2bed.pl` but the BED file is still not compatible so make this file manually in excel or using awk. 

```
split_bam.py -i 1099_S21_L001Aligned.sortedByCoord.out.bam  -r rrna_mod.bed -o output 1> output.txt
```

## Command to run

```
nextflow run main.nf --input 'data/*R{1,2}_001.fastq.gz' --outdir shuaibs_results --ref /home/subudhak/Documents/sara_samples/Iseq_COVID_batch2_and3_corrected_index/RNAgrinder/resources/GRCh38.primary_assembly.genome.fa --gtf /home/subudhak/Documents/sara_samples/Iseq_COVID_batch2_and3_corrected_index/RNAgrinder/resources/gencode.v43.primary_assembly.basic.annotation.gtf --mode PE --rrnaUse ribodetector --index_dir results/00_index/
```
