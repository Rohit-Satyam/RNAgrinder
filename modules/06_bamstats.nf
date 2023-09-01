params.memory = "3g"
params.cpus = 1
params.outdir = "."
params.jobs = 5


process BAMSTATS{
  publishDir "${params.outdir}/06_bamStats", mode: 'copy'
  memory params.memory
  cpus params.cpus
  maxForks params.jobs
        input:
        tuple val(sid), path(bam)
        output:
        path("*")

        script:
  """
  samtools flagstat -@ ${task.cpus} ${bam} > ${sid}.flagstat.tsv
  samtools coverage ${bam} > ${sid}.coverage.txt
  samtools depth -@ ${task.cpus} -s -d 0 -H ${bam} > ${sid}.depth.tsv
  mosdepth -t ${task.cpus} -x ${sid} ${bam.toRealPath()}
  """
}


process RNASEQCINDEX{
  publishDir "${params.outdir}/06_bamStats", mode: 'copy'
  memory params.memory
  cpus params.cpus
  maxForks params.jobs
    output:
    path("gencode.basic.collapsed.gtf")

        script:
  """
  python ${projectDir}/bin/gtexCollapseAnnotation.py  ${params.gtf} gencode.basic.collapsed.gtf
  """
}

process RNASEQC{
  publishDir "${params.outdir}/06_bamStats", mode: 'copy'
  memory params.memory
  cpus params.cpus
  maxForks params.jobs
    input:
      tuple val(sid), path(bam)
      path(gtf)
    output:
      path("*")

        script:
  """
  rnaseqc ${gtf} ${bam} ${params.rnaseqc_ext} .
  """
}

