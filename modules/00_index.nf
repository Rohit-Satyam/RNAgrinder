params.memory = "3g"
params.cpus = 1
params.outdir = "."


process INDEX{
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}", mode: "copy"

    output:
        path("00_index")
    shell:
'''
## Get genome length
        num=$(seqkit fx2tab !{params.ref} -l -n -i | awk -F'\t' '{sum+=$2;} END{print sum;}')
## Compute optimum genomeSAindexNbases parameter
        indexNbases=$(Rscript -e "cat(round(min(14, log2($num)/2 - 1)))")
## No. of contigs/scaffold
        nChr=$(seqkit fx2tab !{params.ref} -l -n -i | wc -l)
## Optimum genomeChrBinNbits
        ChrBin=$(Rscript -e "cat(round(min(18,log2(max(${num}/${nChr},!{params.readLength})))))")
STAR --runThreadN !{task.cpus} --runMode genomeGenerate --genomeDir 00_index \
--genomeFastaFiles !{params.ref} \
--sjdbGTFfile !{params.gtf} --genomeSAindexNbases ${indexNbases} \
--genomeChrBinNbits ${ChrBin} \
!{params.staridx_ext}
    '''
}

