params.memory = "3g"
params.cpus = 1
params.outdir = "."
params.jobs = 1


process RIBODETECTOR{
cpus params.cpus
maxForks params.jobs

publishDir "${params.outdir}/rRNA/ribodetector", mode: 'copy'
        input:
                //tuple val(sid), path(reads)
                tuple val(sid), path(reads)

        output:
                tuple val(sid), path("*.fq.gz")
                path("*.txt")

        shell:
if ("${params.mode}" == "PE")
'''
mean_len=$(bioawk -c fastx '{sum += length($seq)} END {printf "%d", sum / NR}' !{reads[0]})

ribodetector_cpu -t !{task.cpus} -l ${mean_len} -i !{reads[0]} !{reads[1]} -e rrna \
-r !{sid}.rrna_ribodetector.R1.fq !{sid}.rrna_ribodetector.R2.fq --chunk_size 256 \
-o !{sid}.clean_ribodetector.R1.fq !{sid}.clean_ribodetector.R2.fq

raw=$(zcat !{reads[0]}  | wc -l | awk '{print $1/4}')
rrna=$(wc -l !{sid}.rrna_ribodetector.R1.fq | awk '{print $1/4}')
nonrRNA=$(wc -l !{sid}.clean_ribodetector.R1.fq | awk '{print $1/4}')

echo !{sid},$raw,$rrna,$nonrRNA > !{sid}_ribodetector_results.txt
gzip !{sid}.rrna_ribodetector.R1.fq
gzip !{sid}.rrna_ribodetector.R2.fq
gzip !{sid}.clean_ribodetector.R1.fq
gzip !{sid}.clean_ribodetector.R2.fq
'''
else if ("${params.mode}" == "SE")
'''
mean_len=$(bioawk -c fastx '{sum += length($seq)} END {printf "%d", sum / NR}' !{reads})

ribodetector_cpu -t !{task.cpus} -l ${mean_len} -i !{reads} -e rrna \
-r !{sid}.rrna_ribodetector.fq  --chunk_size 256 \
-o !{sid}.clean_ribodetector.fq

raw=$(zcat !{reads}  | wc -l | awk '{print $1/4}')
rrna=$(wc -l !{sid}.rrna_ribodetector.fq | awk '{print $1/4}')
nonrRNA=$(wc -l !{sid}.clean_ribodetector.fq  | awk '{print $1/4}')

echo !{sid},$raw,$rrna,$nonrRNA > !{sid}_ribodetector_results.txt
gzip !{sid}.rrna_ribodetector.fq
gzip !{sid}.clean_ribodetector.fq

'''
}


process SORTMERNA{
cpus params.cpus
maxForks params.jobs
publishDir "$params.outdir/rRNA/sortmerna", mode: 'copy'

        input:
                //tuple val(sid), path(reads)
                tuple val(sid), path(reads)

        output:
                tuple val(sid), path("*.fq.gz")
                path("*.log")
                path("*.txt")

        shell:
if ("${params.mode}" == "PE")
'''
sortmerna !{params.sortmerna_ext} --reads !{reads[0]} --reads !{reads[1]} --fastx --threads !{task.cpus} \
--paired_out --out2 --aligned !{sid}.rrna_sortmerna --other !{sid}.clean_sortmerna --kvdb !{sid}/kvdb \
--idx-dir !{sid}/idx --readb !{sid}/readb

mv !{sid}.rrna_sortmerna_fwd.fq.gz !{sid}.rrna_sortmerna.R1.fq.gz
mv !{sid}.rrna_sortmerna_rev.fq.gz !{sid}.rrna_sortmerna.R2.fq.gz
mv !{sid}.clean_sortmerna_fwd.fq.gz !{sid}.clean_sortmerna.R1.fq.gz
mv !{sid}.clean_sortmerna_rev.fq.gz !{sid}.clean_sortmerna.R2.fq.gz

raw=$(zcat !{reads[0]}  | wc -l | awk '{print $1/4}')
rrna=$(zcat !{sid}.rrna_sortmerna.R1.fq.gz | wc -l | awk '{print $1/4}')
nonrRNA=$(zcat !{sid}.clean_sortmerna.R1.fq.gz | wc -l | awk '{print $1/4}')

echo !{sid},$raw,$rrna,$nonrRNA > !{sid}_sortmerna_results.txt
rm -rf !{sid}/kvdb
'''
else if ("${params.mode}" == "SE")
'''
sortmerna !{params.sortmerna_ext} --reads !{reads[0]} --fastx --threads !{task.cpus} \
--aligned !{sid}.rrna_sortmerna --other !{sid}.clean_sortmerna --kvdb !{sid}/kvdb \
--idx-dir !{sid}/idx --readb !{sid}/readb


mv !{sid}.rrna_sortmerna_fwd.fq.gz !{sid}.rrna_sortmerna.R1.fq.gz
mv !{sid}.rrna_sortmerna_rev.fq.gz !{sid}.rrna_sortmerna.R2.fq.gz
mv !{sid}.clean_sortmerna_fwd.fq.gz !{sid}.clean_sortmerna.R1.fq.gz
mv !{sid}.clean_sortmerna_rev.fq.gz !{sid}.clean_sortmerna.R2.fq.gz

raw=$(zcat !{reads[0]}  | wc -l | awk '{print $1/4}')
rrna=$(zcat !{sid}.rrna_sortmerna.R1.fq.gz | wc -l | awk '{print $1/4}')
nonrRNA=$(zcat !{sid}.clean_sortmerna.R1.fq.gz | wc -l | awk '{print $1/4}')

echo !{sid},$raw,$rrna,$nonrRNA > !{sid}_sortmerna_results.txt
rm -rf !{sid}/kvdb
'''
}

process COMMON_READS{
publishDir "$params.outdir/rRNA/common", mode: 'copy'
cpus params.cpus
maxForks params.jobs

input:
tuple val(sid), path(reads)

output:
tuple val(sid), path("${sid}.R1.fq.gz"), path("${sid}.R2.fq.gz")

shell:
if ("${params.mode}" == "PE")
'''
seqkit common !{reads[0]} !{reads[2]} -s -i -o !{sid}.R1.fq.gz
seqkit common !{reads[1]} !{reads[4]} -s -i -o !{sid}.R2.fq.gz
'''
else if ("${params.mode}" == "SE")
'''
seqkit common !{reads[0]} !{reads[1]} -s -i -o !{sid}.fq.gz
'''
}
