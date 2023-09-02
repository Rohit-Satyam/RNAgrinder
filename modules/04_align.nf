params.memory = "3g"
params.cpus = 1
params.outdir = "."
params.jobs = 1


process STAR{
        publishDir "${params.outdir}/04_alignment/star", mode: 'copy'
        cpus params.cpus
        maxForks params.jobs

input:
        tuple val(sid), path(reads)
        each path(index)

output:
//tuple val(sid), path("*fastq.gz")
tuple val(sid), path("${sid}_Aligned.sortedByCoord.out.bam"), emit: genomic_bam
path("${sid}_Aligned.toTranscriptome.out.bam")
path("*.out"), emit: star_logs
path("*.tab")

shell:
if ("${params.mode}" == "PE")
'''

id=$(zcat !{reads[0]} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//'| awk '{print $1}')
echo ${id}
STAR --genomeDir !{index.toRealPath()} --runThreadN !{task.cpus} \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--chimJunctionOverhangMin 12 \
--chimMultimapNmax 0 \
--chimNonchimScoreDropMin 10 \
--chimOutJunctionFormat 1 \
--chimOutType WithinBAM Junctions SeparateSAMold \
--chimSegmentMin 12 \
--genomeLoad NoSharedMemory \
--outBAMsortingThreadN !{task.cpus} \
--outFileNamePrefix !{sid}_ \
--outFilterMatchNminOverLread 0.3 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.3 \
--outFilterType BySJout \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outReadsUnmapped Fastx \
--outSAMattrRGline ID:${id} SM:!{sid} PL:ILLUMINA \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat \
--readFilesIn !{reads[0]} !{reads[1]} \
--limitBAMsortRAM 1276368436\
--outBAMsortingBinsN 200 \
--sjdbGTFfile !{params.gtf}  \
--sjdbScore 1 2> !{sid}.stderr

## Disabling the output of unmapped reads momentarily
#mv !{sid}_Unmapped.out.mate1 !{sid}_unmapped_R1.fastq
#mv !{sid}_Unmapped.out.mate2 !{sid}_unmapped_R2.fastq
#gzip !{sid}_unmapped_R1.fastq
#gzip !{sid}_unmapped_R2.fastq
'''
else if ("${params.mode}" == "SE")
'''
id=$(zcat !{reads[0]} | head -n 1 | cut -f 3-4 -d":" | sed 's/@//'| awk '{print $1}')
echo ${id}
STAR --genomeDir !{index} --runThreadN !{task.cpus} \
--alignIntronMax 1000000 \
--alignIntronMin 20 \
--alignMatesGapMax 1000000 \
--alignSJDBoverhangMin 1 \
--alignSJoverhangMin 8 \
--chimJunctionOverhangMin 12 \
--chimMultimapNmax 0 \
--chimNonchimScoreDropMin 10 \
--chimOutJunctionFormat 1 \
--chimOutType WithinBAM Junctions SeparateSAMold \
--chimSegmentMin 12 \
--genomeLoad NoSharedMemory \
--outBAMsortingThreadN !{task.cpus} \
--outFileNamePrefix !{sid}_ \
--outFilterMatchNminOverLread 0.3 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 20 \
--outFilterScoreMinOverLread 0.3 \
--outFilterType BySJout \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattrRGline ID:${id} SM:!{sid} PL:ILLUMINA \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat \
--readFilesIn !{reads[0]} \
--sjdbGTFfile !{params.gtf} \
--sjdbScore 1 2> !{sid}.stderr
#mv !{sid}_Unmapped.out.mate !{sid}_unmapped.fastq
#gzip !{sid}_unmapped.fastq
'''
}

process KRAKEN{
	publishDir "${params.outdir}/04_alignment/kraken", mode: 'copy'
	maxForks params.jobs
	cpus params.cpus
	input:
	tuple val(sid), path(reads)
	output:
	file "*"
	script:
        if ("${params.mode}" == "PE")
	"""
	kraken2 --threads ${task.cpus} --confidence 0.1 --paired --db ${params.db}  \
  ${reads[0]} ${reads[1]} --output ${sid} --report-minimizer-data \
  --report ${sid}_kraken_report.txt --gzip-compressed --unclassified-out ${sid}_unclassified#.fasta
  #bracken -r ${params.readlen} -d ${params.db} -i ${sid}_kraken_report.txt -o ${sid}_bracken_report.txt -w ${sid}_kraken_bracken.report
	"""
else if ("${params.mode}" == "SE")
"""
        kraken2 --threads ${task.cpus} --confidence 0.1 --db ${params.db}  \
  ${reads[0]} --output ${sid} --report-minimizer-data \
  --report ${sid}_kraken_report.txt --gzip-compressed --unclassified-out ${sid}_unclassified#.fasta
  #bracken -r ${params.readlen} -d ${params.db} -i ${sid}_kraken_report.txt -o ${sid}_bracken_report.txt -w ${sid}_kraken_bracken.report

"""
}
