#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

if( params.help ) {

log.info """
* -------------------------------------------------
 *  RNAgrinder@KAUST: Analyzing RNASeq Dataset
 * -------------------------------------------------
Usage:
	nextflow run main.nf --input "${params.input}" --outdir ${params.outdir} --ref ${params.ref} \
	--gtf ${params.gtf} --mode ${params.mode} --rrnaUse ${params.rrnaUse} \
	--indexing ${params.indexing}
Input:
	#### Mandatory Arguments ####
	* --gtf: Path to reference GTF file. Default [${params.ref}]
	* --index_dir: Provide path to director where STAR indexes are stored. Use this argument when --indexing is set to false.
	* --input: Path to FastqQ files. Default [${params.input}]
	* --mode: If data is Paired-end pass "PE" else "SE". Default [${params.mode}]
	* --outdir: Path/Name of the output directory. Default [${params.outdir}]
	* --ref: Path to reference fasta file. Default [${params.ref}]
	* --rrnaUse: Choose algorithm to detect rRNA. Possible Values: "ribodetector","sortmerna". Default [${params.rrnaUse}]
	
	Default [${params.index_dir}]

	#### Parameters to pass additional Arguments ####
	* --fastp_ext: Additional arguments to pass to FASTP. Default [${params.fastp_ext}]
	* --fastqc_ext: Additional arguments to pass to FASTQC. Default [${params.fastqc_ext}]
	* --opticalDupPixalDis: Alter the optical duplicate pixel distance according to the Flow-Cell used.
	For details visit https://gatk.broadinstitute.org/hc/en-us/articles/360037224932-MarkDuplicatesSpark. Defult [${params.opticalDupPixalDis}]
	* --rnaseqc_ext: Additional arguments to pass to rnaseqc. Default [${params.rnaseqc_ext}]
	* --sortmerna_ext: Additional arguments to pass to SortMerna. Default [${params.sortmerna_ext}]
	* --staridx_ext: Additional arguments to pass to STAR Indexing Step. Default [${params.staridx_ext}]
	* --staralign_ext:  Additional arguments to pass to STAR ALIGNMENT Step. Default [${params.staralign_ext}]

	#### Parameters to Skip certain Steps ####
	* --skipAlignment: Set this "true" to skip Alignment Step. Default [${params.skipAlignment}]
	* --skipKraken: Set this "true" to skip read classification using Kraken2 and Bracken Step. Default [${params.skipKraken}]
	* --skipRibo:Set this "true" to skip RiboDepletion step (in case of PolyA selection protocol). Default [${params.skipRibo}]
	* --skipTrim: Set this "true" to skip Trimming Step. Default [${params.skipTrim}]
	
	

"""

exit 0
}

include {INDEX} from './modules/00_index'
include {FASTQC} from './modules/01_fastqc'
include {FASTP; POSTTRIMFASTQC; SEQKIT} from './modules/02_fastp'
include {RIBODETECTOR;SORTMERNA; COMMON_READS} from './modules/03_rRNADepletion'
include {STAR;KRAKEN} from './modules/04_align'
include {MARKDUP} from './modules/05_markduplicates'
include {BAMSTATS;RNASEQCINDEX;RNASEQC} from './modules/06_bamstats'
include {MULTIQC as PRETRIM; MULTIQC as POSTTRIM; MULTIQC as SUMMARISEALL;MULTIQC as SORTMERNAQC } from './modules/01_fastqc'

params.help= false
params.input = false
params.outdir= false

params.mode = false


workflow{

if (params.input != false){
	if (params.mode == "PE"){
		Channel.fromFilePairs(params.input, checkIfExists: true )
		.set { input_fastqs }
		} else if (params.mode == "SE") {
			Channel.fromPath(params.input, checkIfExists: true ).map { file -> tuple(file.simpleName, file)}
			.set { input_fastqs }
	}
}
/*
if (params.indexing == true && params.index_dir == "" && params.skipAlignment == false){
	INDEX()
	index=INDEX.out[0]
	} else if (params.indexing == false && params.index_dir == "" && params.skipAlignment == false){
	echo 'Provide a directory containing STAR indexes using --index_dir argument.Otherwise alignments will fail.'
	} else if (params.index_dir != "" && params.skipAlignment == false){
		index = params.index_dir
	}

*/
if (params.index_dir == "" && params.skipAlignment == false){
	INDEX()
	index=INDEX.out[0]
	} else if (params.index_dir != "" && params.skipAlignment == false){
		//index = Channel.fromPath(params.index_dir,  checkIfExists: true, type: 'dir')
		index = Channel.fromPath(params.index_dir)
	}

// Fastqc and Seqkit Summary
	rawfqc_ch=FASTQC(input_fastqs)
	pretrim_input=FASTQC.out.fastqc.collect()
	PRETRIM("01_rawFastQC",pretrim_input,'pre-trimming')
	//SEQKIT()


if (params.skipTrim == false){
        FASTP(input_fastqs)
	FASTP.out[0]
        POSTTRIMFASTQC(FASTP.out[0])
        postrim_input=POSTTRIMFASTQC.out.postfastqc.collect()
        POSTTRIM("02_adapterTrimming",postrim_input,'post-trimming')
		//COMMON_READS()
}

if (params.skipRibo == false && params.skipTrim == true){
	RIBODETECTOR(input_fastqs)} else if (params.skipRibo == false && params.skipTrim == false){
	RIBODETECTOR(FASTP.out[0])}

if (params.skipKraken == false && params.skipTrim == false){
	KRAKEN(FASTP.out[0])} else if (params.skipKraken == false && params.skipTrim == true){
	KRAKEN(input_fastqs)

}

if (params.skipAlignment == false){
	if (params.rrnaUse == "ribodetector" && params.skipRibo == false){
	STAR(RIBODETECTOR.out[0],index)
	MARKDUP(STAR.out.genomic_bam)
	BAMSTATS(MARKDUP.out[0])
	RNASEQCINDEX()
	RNASEQC(MARKDUP.out[0],RNASEQCINDEX.out[0])
	} else if (params.rrnaUse == "sortmerna"){
	SORTMERNA(FASTP.out[0])
	SORTMERNAQC("rRNA/sortmerna",SORTMERNA.out[1].collect(),'post-rRNA-Removal')
	STAR(SORTMERNA.out[0],index)
	MARKDUP(STAR.out.genomic_bam)
	BAMSTATS(MARKDUP.out[0])
	RNASEQCINDEX()
	RNASEQC(MARKDUP.out[0],RNASEQCINDEX.out[0])
	} else {
	STAR(input_fastqs,index)
	MARKDUP(STAR.out.genomic_bam)
	BAMSTATS(MARKDUP.out[0])
	RNASEQCINDEX()
	RNASEQC(MARKDUP.out[0],RNASEQCINDEX.out[0])
	}
}
if (params.skipTrim == true && params.skipAlignment == true ){
	all_combine=SORTMERNA.out[1]
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')
	} else if (params.skipTrim == true && params.skipAlignment == false ){
	all_combine=BAMSTATS.out[0].mix(MARKDUP.out.dedupmtx,STAR.out.star_logs,RNASEQC.out[0])
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')
	} else if (params.skipTrim == false && params.skipAlignment == true ){
	all_combine=FASTP.out.fastp_logs
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')
	} else {
	all_combine=BAMSTATS.out[0].mix(MARKDUP.out.dedupmtx,STAR.out.star_logs,RNASEQC.out[0],FASTP.out.fastp_logs)
	SUMMARISEALL("Summary_Reports",all_combine.collect(),'summary_report')
	}


}
