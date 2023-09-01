params.memory = "3g"
params.cpus = 1
params.outdir = "."
params.jobs = 5


process MARKDUP{
  publishDir "${params.outdir}/05_markDuplicates", pattern: "*.dedup.{bam,bai}", mode: 'copy'
  publishDir "${params.outdir}/05_markDuplicates", pattern: "*.txt", mode: 'copy'

  memory params.memory
  cpus params.cpus
  maxForks params.jobs
  input:
  tuple val(sid), path(bam)
  output:
  tuple val(sid), path("${sid}.dedup.bam")
  path("${sid}_markdup_metrics.txt"), emit: dedupmtx
  path("${sid}.dedup.bai")
  
  script:
  """
  picard MarkDuplicates I=${bam} O=${sid}.dedup.bam M=${sid}_markdup_metrics.txt \
  CREATE_INDEX=true READ_NAME_REGEX=nulll ASSUME_SORT_ORDER=coordinate
  #gatk MarkDuplicatesSpark -I ${bam} -O ${sid}.dedup.bam -M ${sid}_markdup_metrics.txt --tmp-dir . -OBI
  """
}
