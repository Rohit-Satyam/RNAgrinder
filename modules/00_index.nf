params.memory = "3g"
params.cpus = 1
params.outdir = "."


process INDEX{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}", mode: "copy"

    output:
        path("00_index")
    """
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir 00_index \
--genomeFastaFiles ${params.ref} \
--sjdbGTFfile ${params.gtf} ${params.staridx_ext}
    """
}

