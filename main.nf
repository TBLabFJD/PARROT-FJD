#!/usr/bin/env nextflow
import java.util.regex.Matcher

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl=2



/*
 * Print the default parameters
 */ 

log.info """\
FJD-PIPELINE    -    NEXTFLOW  
=============================
samples      : $params.samples
"""


/* 
 * Import modules 
 */
include { BS_CHECK } from './modules/execution_modules'
include { LOCAL_CHECK } from './modules/execution_modules'
include { BS_COPY } from './modules/execution_modules'
include { FASTQ_CONCATENATION } from './modules/execution_modules'
include { FASTP } from './modules/execution_modules' //YBQ: Fastp for splut fastqs
include { BWA } from './modules/execution_modules'
include { FASTQTOSAM } from './modules/execution_modules'
include { MERGEBAMALIGNMENT } from './modules/execution_modules'
include { MARKDUPLICATESSPARK } from './modules/execution_modules'
include { SORTSAM } from './modules/execution_modules'
include { SETTAGS } from './modules/execution_modules'
include { BASERECALIBRATOR } from './modules/execution_modules'
include { APPLYBQSR } from './modules/execution_modules'
include { MERGEBAM } from './modules/execution_modules' //YBQ: Merge bams when mapping run in parallel
include { LOCALBAM as LOCALBAM } from './modules/execution_modules'
include { LOCALBAM as LOCALBAM_CNV } from './modules/execution_modules'

//nuevo modulo GUR: SPLIT_BAM y MERGE_SPLIT_VCF
include { SPLIT_BAM } from './modules/execution_modules'

////GUR: nuevos modulos merge vcfs por cromosomas de cada programa
include { MERGE_SPLIT_VCF as MERGE_SPLIT_VCF_GATK} from './modules/execution_modules'
include { MERGE_SPLIT_VCF as MERGE_SPLIT_VCF_DEEPVARIANT } from './modules/execution_modules'
include { MERGE_SPLIT_VCF as MERGE_SPLIT_VCF_DRAGEN } from './modules/execution_modules'

include { PARALLEL_HAPLOTYPECALLER } from './modules/execution_modules'

include { HAPLOTYPECALLER } from './modules/execution_modules'
include { SELECT_SNV } from './modules/execution_modules'
include { SELECT_INDEL } from './modules/execution_modules'
include { SELECT_MIX } from './modules/execution_modules'
include { FILTRATION_SNV } from './modules/execution_modules'
include { FILTRATION_INDEL } from './modules/execution_modules'
include { FILTRATION_MIX } from './modules/execution_modules'
include { MERGE_VCF } from './modules/execution_modules'
include { DEEPVARIANT } from './modules/execution_modules'
// nuevo módulo PARALLEL_DRAGEN: PARALELIZA Y ADEMAS tiene los dos pasos de DRAGEN en 1: STR_MODEL + HAPLOTYPE_CALLER DRAGEN 
include { PARALLEL_HAPLOTYPECALLER_DRAGEN } from './modules/execution_modules'

include { STR_MODEL_DRAGEN } from './modules/execution_modules'
include { HAPLOTYPECALLER_DRAGEN } from './modules/execution_modules'
include { FILTRATION_DRAGEN } from './modules/execution_modules'
include { FILTER_VCF as FILTER_VCF_GATK } from './modules/execution_modules'
include { FILTER_VCF as FILTER_VCF_DEEPVARIANT } from './modules/execution_modules'
include { FILTER_VCF as FILTER_VCF_DRAGEN } from './modules/execution_modules'
//GUR: mis nuevos modulos
include { FINAL_VCF as FINAL_GATK } from './modules/execution_modules'
include { FINAL_VCF as FINAL_DRAGEN } from './modules/execution_modules'
include { FINAL_VCF as FINAL_DEEPVARIANT } from './modules/execution_modules'
/// GUR
include { BAM2CRAM } from './modules/execution_modules' 
include { CRAM2BAM } from './modules/execution_modules' 


///
include { MERGE_VCF_CALLERS } from './modules/execution_modules'

include { GVCF_HAPLOTYPECALLER } from './modules/execution_modules'
include { COMBINE_GVCF } from './modules/execution_modules'
include { GENOTYPE_GVCF } from './modules/execution_modules'
include { CALCULATE_GENOTYPE_POSTERIORS } from './modules/execution_modules'
include { VARIANT_FILTRATION } from './modules/execution_modules'
include { VARIANT_ANNOTATOR } from './modules/execution_modules'

include { LOCALVCF } from './modules/execution_modules'
include { FORMAT2INFO } from './modules/execution_modules'
include { AUTOMAP } from './modules/execution_modules'
include { VEP } from './modules/execution_modules'
include { PVM } from './modules/execution_modules'

include { MOSDEPTH } from './modules/execution_modules'
//include { MOSDEPTH_JOIN as MOSDEPTH_JOIN_SNV } from './modules/execution_modules'
//include { MOSDEPTH_JOIN as MOSDEPTH_JOIN_CNV } from './modules/execution_modules'
include { MOSDEPTH_PLOT } from './modules/execution_modules'
include { MOSDEPTH_COV } from './modules/execution_modules'
include { GENOMECOV } from './modules/execution_modules'
include { SAMTOOLS_FLAGSTAT } from './modules/execution_modules'
include { READ_LENGTH_STATS } from './modules/execution_modules'
include { SEQUENCING_QUALITY_SCORES } from './modules/execution_modules'
include { SEQUENCING_CG_AT_CONTENT } from './modules/execution_modules'
include { NREADS_NONDUP_UNIQ } from './modules/execution_modules'
include { QC_SUMMARY } from './modules/execution_modules'
include { RUN_QC_CAT } from './modules/execution_modules'

include { BEDPROCCESING } from './modules/execution_modules'
include { EXOMEDEPTH } from './modules/execution_modules'
include { CONVADING } from './modules/execution_modules'
include { PANELCNMOPS } from './modules/execution_modules'
include { CNV_RESULT_MIXER } from './modules/execution_modules'
include { ANNOTSV } from './modules/execution_modules'
include { PAM } from './modules/execution_modules'

include { MANTA } from './modules/execution_modules'
include { ANNOTSV_VCF } from './modules/execution_modules'

include { PREPROCESSINTERVALS } from './modules/execution_modules'
include { COLLECTREADCOUNTS } from './modules/execution_modules'
include { ANNOTATEINTERVALS } from './modules/execution_modules'
include { FILTERINTERVALS } from './modules/execution_modules'
include { INTERVALLISTTOOLS } from './modules/execution_modules'
include { DETERMINEGERMLINECONTIGPLOIDY } from './modules/execution_modules'
include { GERMLINECNVCALLER } from './modules/execution_modules'
include { POSTPROCESSGERMLINECNVCALLS } from './modules/execution_modules'
include { VCF2BED } from './modules/execution_modules'



// def final_vcf  = Channel.fromPath(params.final_vcf)
// def omim       = params.omim       ? Channel.fromPath(params.omim)       : Channel.empty()
// def genefilter = params.genefilter ? Channel.fromPath(params.genefilter) : Channel.empty()
// def glowgenes  = params.glowgenes  ? Channel.fromPath(params.glowgenes)  : Channel.empty()

// println "final_vcf: ${final_vcf}"
// // println "${omim}"

// checkPathParamList = [
//     params.dbNSFP_gene, params.omim,
//     params.genefilter, params.glowgenes
// ]
// // for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
// if (params.dbNSFP_gene) { ch_input = Channel.fromPath(params.dbNSFP_gene) } else { exit 1, 'dbNSFP_gene samplesheet not specified!' }
// if (genefilter) { println "HOLAAAAAAAA" } else { exit 1, 'genefilter samplesheet not specified!' }
// println "dbNSFP_gene: ${params.dbNSFP_gene}"
// println "omim: ${omim}"
// println omim.view()
// println "genefilter: ${genefilter}"
// println "glowgenes: ${params.glowgenes}"

// if (params.genefilter) { println "ADIOS" } else { exit 1, 'genefilter samplesheet not specified!' }

// println Channel.fromPath(params.dbNSFP_gene).type()
// println file(params.omim).type()

// exit 1

// aa = params.analysis.toUpperCase().split(",")

// println aa


// if(aa.contains("M"))  {println "HOLA"} else {println "ADIOS"}
// if(params.analysis.toUpperCase().contains("M"))  {println "1"} else {println "2"}

// channel.value(params.analysis).view()





  
workflow CHECK_PARAMS {

	take:


	main:
		println "Checking..."

		//Check the existance of reference genome
		if ( !params.reference_fasta ) {exit 1, "Error: Missing reference fasta file definition.\n"} 
		else {file(params.reference_fasta, type: "file", checkIfExists: true)}

		if ( !params.reference_index ) {exit 1, "Error: Missing index reference fasta file definition.\n"} 
		else {file(params.reference_index, type: "file", checkIfExists: true)}

		if ( !params.reference_dict  ) {exit 1, "Error: Missing dictionary reference fasta file definition.\n"} 
		else {file(params.reference_dict, type: "file", checkIfExists: true)}

		if ( !params.reference_gzi   ) {exit 1, "Error: Missing gzi reference fasta file definition.\n"} 
		else {file(params.reference_gzi, type: "file", checkIfExists: true)}
		println "Referencie file check"



		// Check the existance of the input file
		if ( !params.input ) {exit 1, "Error: Missing input file definition.\n"}
		if ( params.input && !params.analysis.toUpperCase().contains("D") ) { file(params.input, type: "dir", checkIfExists: true) }
		println "Input folder check"
		


		// Define the run name
		if (params.runname) { runname = params.runname }
		else if ( params.analysis.toUpperCase().contains("D") )  { runname = params.input }
		else { runname = new Date().format("yyyy-MM-dd_HH-mm") }
		println "Run name: $runname" 



		// Check the existance of the output directory
		if (params.output) { 
			file(params.output, type: "dir").mkdir()
			outputdir = file(params.output, type: "dir", checkIfExists: true) 
		} else {
			exit 1, "Error: Missing output directory definition.\n"
		}
		println "outputdir check: $outputdir"



		// Check the existance of the sample file
		if ( params.samples ) { file(params.samples, type: "file", checkIfExists: true) }
		println "Sample file check"



		// Check the existance of the bed file
		// cnv_analysis_check = params.analysis.toUpperCase() =~ /[CX]/
		// assert cnv_analysis_check instanceof Matcher
		// if ( !params.bed && ( cnv_analysis_check || params.intervals ) && params.capture != "G") {
		if ( !params.bed && ( (params.analysis.toUpperCase() =~ /[CX]/) || params.intervals ) && params.capture != "G") {
			exit 1, "Error: No bedfile provided and CNV calling or '--intervals' was provided for Panels or WES analisis.\n"
		} else if ( params.bed ) {
			bedfile = file(params.bed, type: "file", checkIfExists: true)
			bedfilecontent = bedfile.readLines()[0].split("\t").length
			if ( bedfilecontent != 4 ) {exit 1, "Error: Bedfile must have 4 columns.\n"}
		}
		println "Bed file check"



		// Check the existance of the ped (pedigree) file
		if ( params.ped &&  params.analysis.toUpperCase().contains("G") ) { pedfile = file(params.ped, type: "file", checkIfExists: true) }
		if ( params.ped && !params.analysis.toUpperCase().contains("G") ) { println "WARN: Ped (pedigree) file provided, but not used" }
		println "Ped file check"



		// Check the type of the analysis
		// m = params.analysis.toUpperCase() =~ /[DMQSGCXAN]/
		// assert m instanceof Matcher
		// if ( !m ) {exit 1, "Error: Cannot recognice the any of the specified analisis analysis.\nThe available analysis are: D (Download from BaseSpace), M (Mapping), S (SNV individual), G (SNV GVCF),\nA (SNV annotation), C (CNV calling), N (CNV annotation), X (chrX CNV calling)\n"}
		if ( !(params.analysis.toUpperCase() =~ /[DMQSGCXAN]/) ) {exit 1, "Error: Cannot recognice the any of the specified analisis analysis.\nThe available analysis are: D (Download from BaseSpace), M (Mapping), S (SNV individual), G (SNV GVCF),\nA (SNV annotation), C (CNV calling), N (CNV annotation), X (chrX CNV calling)\n"}
		println "Analysis type check"



		// BaseSpace must start with mapping
		/*if ( !params.analysis.toUpperCase().contains("M") && params.analysis.toUpperCase().contains("D") ) {exit 1, "Error: If basespace parameter is specify, the mapping (M) analysis must be specify.\n"}
		println "Incompatibility check"*/



		if ( params.analysis.toUpperCase().contains("D") ) {

			BS_CHECK(
				params.input,
				params.baseuser,
				params.samples,
				params.analysis.toUpperCase() )
			controlsamples  = BS_CHECK.out.controlsamples
			samples2analyce = BS_CHECK.out.samples2analyce
			datasets = BS_CHECK.out.datasets


		} else {

			LOCAL_CHECK(
				params.input,
				params.samples,
				params.analysis.toUpperCase() )

			controlsamples  = LOCAL_CHECK.out.controlsamples
			samples2analyce = LOCAL_CHECK.out.samples2analyce
			datasets = []
		}

		

	emit:
		controlsamples_file = controlsamples
		samples2analyce_file = samples2analyce
		controlsamples  = controlsamples.splitCsv()
		samples2analyce = samples2analyce.splitCsv()
		runname = runname
		datasets = datasets
}



workflow DOWNLOAD {
	take:
		samplename
		datasets

	main:

		if ( params.analysis.toUpperCase().contains("D") ) {

			BS_COPY (
				params.input,
				samplename,
				params.baseuser,
				datasets )

			fastq = BS_COPY.out.fastq

		} else {

			FASTQ_CONCATENATION (
				params.input,
				samplename )

			fastq = FASTQ_CONCATENATION.out.fastq
		}

	emit:
		fastq = fastq
}







workflow MAPPING {
	take:
		fastq

	main:

		BWA (
			fastq,
			params.reference_fasta,
			params.bwa_amb,
			params.bwa_ann,
			params.bwa_pac,
			params.bwa_bwt,
			params.bwa_sa )

		FASTQTOSAM (
			fastq,
			params.scratch )

		
		MERGEBAMALIGNMENT (
			BWA.out.mapped_bam.join(FASTQTOSAM.out.unmapped_bam),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )
		
		if ( !params.ignoreduplicates ) {
			
			MARKDUPLICATESSPARK (
				MERGEBAMALIGNMENT.out.merged_bam,
				params.scratch )

			sorted_bam = MARKDUPLICATESSPARK.out.deduppedsorted_bam

		} else {

			SORTSAM (
				MERGEBAMALIGNMENT.out.merged_bam,
				params.scratch )

			sorted_bam = SORTSAM.out.sorted_bam
		}


		SETTAGS (
			sorted_bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )


		BASERECALIBRATOR (
			SETTAGS.out.tagged_bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.g1000_knownsites,
			params.g1000_knownsites_idx,
			params.mills_knownsites,
			params.mills_knownsites_idx,
			params.dbsnp_knownsites,
			params.dbsnp_knownsites_idx,
			params.scratch )


		APPLYBQSR (
			SETTAGS.out.tagged_bam.join(BASERECALIBRATOR.out.bqsr_table),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch,
			params.assembly )


	emit:
		bam = APPLYBQSR.out.bam

}






workflow SNVCALLING {
	take:
		bam
	main:

		HAPLOTYPECALLER (
			bam,
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_SNV (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			HAPLOTYPECALLER.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTER_VCF_GATK (
			MERGE_VCF.out.vcf,
			params.assembly,
			"gatk" )

		


		DEEPVARIANT (
			bam,
			params.bed,
			params.intervals,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.capture,
			params.scratch )

		FILTER_VCF_DEEPVARIANT (
			DEEPVARIANT.out.vcf,
			params.assembly,
			"deepvariant" )




		STR_MODEL_DRAGEN(
			bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.reference_str,
			params.scratch )

		HAPLOTYPECALLER_DRAGEN(
			bam.join(STR_MODEL_DRAGEN.out.strmodel),
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_DRAGEN(
			HAPLOTYPECALLER_DRAGEN.out.vcf,
			params.scratch )

		FILTER_VCF_DRAGEN(
			FILTRATION_DRAGEN.out.vcf,
			params.assembly,
			"dragen" )

		



		MERGE_VCF_CALLERS(
			FILTER_VCF_GATK.out.vcf.join(FILTER_VCF_DEEPVARIANT.out.vcf).join(FILTER_VCF_DRAGEN.out.vcf),
			params.assembly,
			params.reference_fasta,
			projectDir )




	emit:
		finalvcf = MERGE_VCF_CALLERS.out.vcf
}


////////////////////////////////////////////////////// START: nuevos workflows GUR 9 febrero 2023 ///////////////////////////////

workflow PARALLEL_GATKCALLING  {
	take:
		bam
	main:

		PARALLEL_HAPLOTYPECALLER (
			bam,
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch 
		)
		MERGE_SPLIT_VCF_GATK (
			PARALLEL_HAPLOTYPECALLER.out.vcf.groupTuple(),
			params.reference_fasta,
			params.scratch,
			"gatk"
		)
		
	emit:
		finalvcf = MERGE_SPLIT_VCF_GATK.out.vcf
}

workflow PROCESS_GATKCALLING {
	take:
		vcf //(y sample??) -> mirar donde take bam porque take bam de normal era una tupla que traia ya el sample
	main:
		SELECT_SNV (
			vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTER_VCF_GATK (
			MERGE_VCF.out.vcf,
			params.assembly,
			"gatk" )

		// nuevo proceso GUR: este proceso es el FINAL_VCF que lo llamamos con FINAL_GATK DRAGEN Y DEEP VARIANT INDEPENDIENTE///
		FINAL_GATK(
			FILTER_VCF_GATK.out.vcf,
			params.assembly,
			"gatk" )

	emit:
		finalvcf = FINAL_GATK.out.vcf
}

workflow PARALLEL_DRAGENCALLING {
	take:
		bam
	main:
		PARALLEL_HAPLOTYPECALLER_DRAGEN(
			bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.reference_str,
			params.scratch,			
			params.bed,
			params.intervals,
			params.padding )
		
		MERGE_SPLIT_VCF_DRAGEN (
			PARALLEL_HAPLOTYPECALLER_DRAGEN.out.vcf.groupTuple(),
			params.reference_fasta,
			params.scratch,
			"dragen"
		)

	emit:
		finalvcf = MERGE_SPLIT_VCF_DRAGEN.out.vcf
}

workflow PROCESS_DRAGENCALLING {
	take:
		vcf
	main:
		FILTRATION_DRAGEN(
			vcf,
			params.scratch )

		FILTER_VCF_DRAGEN(
			FILTRATION_DRAGEN.out.vcf,
			params.assembly,
			"dragen" )

		// nuevo proceso GUR: este proceso es el FINAL_VCF que lo llamamos con FINAL_GATK DRAGEN Y DEEP VARIANT INDEPENDIENTE///
		FINAL_DRAGEN(
			FILTER_VCF_DRAGEN.out.vcf,
			params.assembly,
			"dragen" )
			
	emit:
		finalvcf = FINAL_DRAGEN.out.vcf
}



//nuevo workflow YO -> DEEPVARIAN(T)
workflow DEEPVARIANTCALLING {
	take:
		bam
	main:

		SPLIT_BAM (
			bam )

		DEEPVARIANT (
			SPLIT_BAM.out.transpose(),
			params.bed,
			params.intervals,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.capture,
			params.scratch )
		
		MERGE_SPLIT_VCF_DEEPVARIANT (
			DEEPVARIANT.out.vcf.groupTuple(),
			params.reference_fasta,
			params.scratch,
			"deepvariant"
		)

		FILTER_VCF_DEEPVARIANT (
			MERGE_SPLIT_VCF_DEEPVARIANT.out.vcf,
			params.assembly,
			"deepvariant" )

		// nuevo proceso GUR: este proceso es el FINAL_VCF que lo llamamos con FINAL_GATK DRAGEN Y DEEP VARIANT INDEPENDIENTE///
		FINAL_DEEPVARIANT(
			FILTER_VCF_DEEPVARIANT.out.vcf,
			params.assembly,
			"deepvariant" )
			
	emit:
		finalvcf = FINAL_DEEPVARIANT.out.vcf
}



// workflow GATK calling parallel TODO: ya es antiguo
/*workflow GATKCALLING {
	take:
		bam
	main:

		SPLIT_BAM (
			bam )

		PARALLEL_HAPLOTYPECALLER (
			SPLIT_BAM.out.transpose(),
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch 
		)
		MERGE_SPLIT_VCF_GATK (
			HAPLOTYPECALLER.out.vcf.groupTuple(),
			params.reference_fasta,
			params.scratch,
			"gatk"
		)
		SELECT_SNV (
			MERGE_SPLIT_VCF_GATK.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			MERGE_SPLIT_VCF_GATK.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			MERGE_SPLIT_VCF_GATK.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTER_VCF_GATK (
			MERGE_VCF.out.vcf,
			params.assembly,
			"gatk" )
		
		// nuevo proceso GUR: este proceso es el FINAL_VCF que lo llamamos con FINAL_GATK DRAGEN Y DEEP VARIANT INDEPENDIENTE///
		FINAL_GATK(
			FILTER_VCF_GATK.out.vcf,
			params.assembly,
			"gatk" )

	emit:
		finalvcf = FINAL_GATK.out.vcf
}*/

//nuevo workflow YO -> DRAG(E)N
/* workflow DRAGENCALLING {
	take:
		bam
	main:
		SPLIT_BAM (
			bam )

		PARALLEL_HAPLOTYPECALLER_DRAGEN(
			SPLIT_BAM.out.transpose(),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.reference_str,
			params.scratch,			
			params.bed,
			params.intervals,
			params.padding )
		
		MERGE_SPLIT_VCF_DRAGEN (
			PARALLEL_HAPLOTYPECALLER_DRAGEN.out.vcf.groupTuple(),
			params.reference_fasta,
			params.scratch,
			"dragen"
		)

		FILTRATION_DRAGEN(
			MERGE_SPLIT_VCF_DRAGEN.out.vcf,
			params.scratch )

		FILTER_VCF_DRAGEN(
			FILTRATION_DRAGEN.out.vcf,
			params.assembly,
			"dragen" )

		// nuevo proceso GUR: este proceso es el FINAL_VCF que lo llamamos con FINAL_GATK DRAGEN Y DEEP VARIANT INDEPENDIENTE///
		FINAL_DRAGEN(
			FILTER_VCF_DRAGEN.out.vcf,
			params.assembly,
			"dragen" )
			
	emit:
		finalvcf = FINAL_DRAGEN.out.vcf
} */

////////////////////////////////////////////////////// END: nuevos workflows GUR 9 febrero 2023 ///////////////////////////////


workflow ANNOTATION {

	take:
		final_vcf

	main:

		FORMAT2INFO( 
			final_vcf )


		AUTOMAP(
			final_vcf,
			params.assembly,
			projectDir )


		VEP(
			params.dbscSNV,
			params.dbscSNV_tbi,
			params.loFtool,
			params.exACpLI,
			params.dbNSFP,
			params.dbNSFP_tbi,
			params.maxEntScan,
			params.cADD_INDELS,
			params.cADD_INDELS_tbi,
			params.cADD_SNVS,
			params.cADD_SNVS_tbi,
			params.kaviar,
			params.kaviar_tbi,
			params.cCRS_DB,
			params.cCRS_DB_tbi,
			params.dENOVO_DB,
			params.dENOVO_DB_tbi,
			params.cLINVAR,
			params.cLINVAR_tbi,
			params.gNOMADg,
			params.gNOMADg_tbi,
			params.gNOMADe,
			params.gNOMADe_tbi,
			params.gNOMADg_cov,
			params.gNOMADg_cov_tbi,
			params.gNOMADe_cov,
			params.gNOMADe_cov_tbi,
			params.cSVS,
			params.cSVS_tbi,
			params.mutScore,
			params.mutScore_tbi,
			params.mAF_FJD_COHORT,
			params.mAF_FJD_COHORT_tbi,
			params.spliceAI_SNV,
			params.spliceAI_SNV_tbi,
			params.spliceAI_INDEL,
			params.spliceAI_INDEL_tbi,
			params.vep_cache,
			params.vep_plugins,
			params.vep_fasta,
			params.vep_fai,
			params.vep_gzi,
			params.vep_assembly,
			final_vcf.join(FORMAT2INFO.out.sample_info),
			params.assembly )



		PVM(
			VEP.out.vep_tsv.join(AUTOMAP.out.roh_automap),
			params.dbNSFP_gene,
			params.omim,
			params.regiondict,
			params.maf,
			params.genelist,
			params.glowgenes,
			params.assembly,
			projectDir )


	// emit:


}







workflow CNVCALLING {
	take:
		all_bam
		all_bai
		runname
		samples2analyce
		bam

	main:

		// Bed file transformation
		BEDPROCCESING(
			params.bed,
			params.min_target,
			params.window,
			params.cnv_chrx,
			params.fai_convading,
			projectDir )

		// BEDPROCCESING2(
		// 	params.bed,
		// 	params.min_target,
		// 	params.window,
		// 	params.cnv_chrx,
		// 	params.fai_convading,
		// 	// projectDir,
		// 	params.reference_fasta,
		// 	params.reference_index,
		// 	params.reference_dict,
		// 	params.reference_gzi )


		// CNV calling
		EXOMEDEPTH(
			all_bam,
			all_bai,
			BEDPROCCESING.out.bed,
			runname,
			projectDir )

		CONVADING(
			all_bam,
			all_bai,
			BEDPROCCESING.out.bed,
			runname,
			params.fai_convading,
			projectDir )
		
		PANELCNMOPS(
			all_bam,
			all_bai,
			BEDPROCCESING.out.bed,
			runname,
			projectDir )


		// // GATK CNV calling
		// PREPROCESSINTERVALS (
		// 	params.reference_fasta,
		// 	params.reference_index,
		// 	params.reference_dict,
		// 	params.reference_gzi,
		// 	params.bed,
		// 	runname )


		// COLLECTREADCOUNTS (
		// 	bam,
		// 	params.reference_fasta,
		// 	params.reference_index,
		// 	params.reference_dict,
		// 	params.reference_gzi,
		// 	PREPROCESSINTERVALS.out.interval_list )


		// ANNOTATEINTERVALS (
		// 	params.reference_fasta,
		// 	params.reference_index,
		// 	params.reference_dict,
		// 	params.reference_gzi,
		// 	PREPROCESSINTERVALS.out.interval_list )

		
		// FILTERINTERVALS (
		// 	PREPROCESSINTERVALS.out.interval_list,
		// 	ANNOTATEINTERVALS.out.annotated_intervals,
		// 	COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList() )

		
		// INTERVALLISTTOOLS(
		// 	FILTERINTERVALS.out.filter_intervals )


		// DETERMINEGERMLINECONTIGPLOIDY(
		// 	FILTERINTERVALS.out.filter_intervals,
		// 	COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList(),
		// 	params.contig_ploidy_priors )


		// GERMLINECNVCALLER(
		// 	DETERMINEGERMLINECONTIGPLOIDY.out.ploidy_calls,
		// 	ANNOTATEINTERVALS.out.annotated_intervals,
		// 	COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList(),
		// 	INTERVALLISTTOOLS.out.scatterout.flatten() )


		// POSTPROCESSGERMLINECNVCALLS(
		// 	DETERMINEGERMLINECONTIGPLOIDY.out.ploidy_calls,
		// 	DETERMINEGERMLINECONTIGPLOIDY.out.sample_index_list.splitCsv(),
		// 	GERMLINECNVCALLER.out.gcnv_dir.collect().flatten().filter( ~/.*calls$/ ).toList(),
		// 	GERMLINECNVCALLER.out.gcnv_dir.collect().flatten().filter( ~/.*model$/ ).toList(),
		// 	params.reference_dict )


		// VCF2BED(
		// 	POSTPROCESSGERMLINECNVCALLS.out.cnv_out.collect().flatten().filter( ~/genotyped-intervals.*/ ).toList(),
		// 	runname )








		// // // CNV merge
		// CNV_RESULT_MIXER(
		// 	EXOMEDEPTH.out.cnvs.join(CONVADING.out.cnvs).join(PANELCNMOPS.out.cnvs).join(VCF2BED.out.cnvs),
		// 	samples2analyce,
		// 	projectDir )
		
		CNV_RESULT_MIXER(
			EXOMEDEPTH.out.cnvs.join(CONVADING.out.cnvs).join(PANELCNMOPS.out.cnvs),
			samples2analyce,
			projectDir )


		// CNV ANNOTATION
		ANNOTSV ( 
			CNV_RESULT_MIXER.out.merged_bed,
			params.genelist,
			params.annotsv_assembly,
			params.annotsv_path )

		PAM (
			ANNOTSV.out.annotated_cnv,
			CNV_RESULT_MIXER.out.colnames,
			params.genelist,
			params.glowgenes,
			projectDir )

	// emit:

}









workflow COMBINEDSNVCALLING {
	take:
		bam
		runname
	main:

		GVCF_HAPLOTYPECALLER (
			bam,
			params.bed,
			params.intervals,
			params.padding,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		vcf = GVCF_HAPLOTYPECALLER.out.gvcf.collect().flatten().filter( ~/.*g.vcf$/ ).toList()
		idx = GVCF_HAPLOTYPECALLER.out.gvcf.collect().flatten().filter( ~/.*g.vcf.idx$/ ).toList()

		COMBINE_GVCF(
			vcf,
			idx,
			runname,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		GENOTYPE_GVCF(
			COMBINE_GVCF.out.gvcf,
			params.ped,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )




		SELECT_SNV (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_INDEL (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		SELECT_MIX (
			GENOTYPE_GVCF.out.vcf,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )

		FILTRATION_SNV (
			SELECT_SNV.out.vcf,
			params.scratch )

		FILTRATION_INDEL (
			SELECT_INDEL.out.vcf,
			params.scratch )

		FILTRATION_MIX (
			SELECT_MIX.out.vcf,
			params.scratch )

		MERGE_VCF (
			FILTRATION_SNV.out.vcf.join(FILTRATION_INDEL.out.vcf).join(FILTRATION_MIX.out.vcf),
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.scratch )


		if (params.ped){

			CALCULATE_GENOTYPE_POSTERIORS(
				MERGE_VCF.out.vcf,
				params.ped,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			VARIANT_FILTRATION(
				CALCULATE_GENOTYPE_POSTERIORS.out.vcf,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			VARIANT_ANNOTATOR(
				VARIANT_FILTRATION.out.vcf,
				params.ped,
				params.reference_fasta,
				params.reference_index,
				params.reference_dict,
				params.reference_gzi,
				params.scratch )

			vcf2filt = VARIANT_ANNOTATOR.out.vcf
			
		} else {
				
			vcf2filt = MERGE_VCF.out.vcf
		}


		FILTER_VCF_GATK (
			vcf2filt,
			params.assembly,
			"gatk" )

	emit:
		finalvcf = FILTER_VCF.out.vcf
}










workflow QUALITYCHECK {
	take:
		bam
		runname
	
	main:

		MOSDEPTH(
			bam,
			params.bed )

		MOSDEPTH_PLOT(
			MOSDEPTH.out.mosdepth.collect(),
			projectDir )

		MOSDEPTH_COV(
			bam,
			params.bed,
			params.padding )



		// // SNV quality check
		// // if ( params.analysis.toUpperCase().contains("S") ) {
		// if ( params.analysis.toUpperCase().contains("Q") ) {

		// 	MOSDEPTH_JOIN_SNV(
		// 		MOSDEPTH.out.mosdepth.collect(),
		// 		params.qsnvthreshold,
		// 		"snv",
		// 		runname )
		// }



		// // CNV quality check
		// // if ( params.analysis.toUpperCase().contains("C") ) {
		// if ( params.analysis.toUpperCase().contains("Q") ) {
			
		// 	MOSDEPTH_JOIN_CNV(
		// 		MOSDEPTH.out.mosdepth.collect(),
		// 		params.qcnvthreshold,
		// 		"cnv",
		// 		runname )
		// }




		GENOMECOV(
			bam,
			params.bed )

		SAMTOOLS_FLAGSTAT(
			bam )

		READ_LENGTH_STATS(
			bam )

		SEQUENCING_QUALITY_SCORES(
			bam )

		SEQUENCING_CG_AT_CONTENT(
			bam )

		NREADS_NONDUP_UNIQ(
			bam )

		QC_SUMMARY(
			GENOMECOV.out.genomecov.join(SAMTOOLS_FLAGSTAT.out.flagstat).join(READ_LENGTH_STATS.out.read_length_stats).join(SEQUENCING_QUALITY_SCORES.out.quality).join(SEQUENCING_CG_AT_CONTENT.out.cg_at).join(NREADS_NONDUP_UNIQ.out.nreads_nondup_uniq),
			params.bed,
			projectDir	)


		// library_stats = QC_SUMMARY.out.quality_summary.collect().flatten().filter( ~/.*library.stats.txt$/ ).toList()
		quality_summary = QC_SUMMARY.out.quality_summary.collect().flatten().filter( ~/.*quality.summary.txt$/ ).toList()

		RUN_QC_CAT(
			quality_summary,
			runname )
	// emit:
}




workflow CNVCALLING_WGS {
	take:
		bam
		bai
		runname
	
	main:

		MANTA (
			bam,
			bai,
			params.reference_fasta,
			params.reference_index,
			params.reference_gzi,
			runname )

		ANNOTSV_VCF (
			MANTA.out.diploidsv,
			params.annotsv_assembly,
			params.annotsv_path )

	// emit:

}





workflow CNVCALLING_WES {
	take:
		bam
		runname
	
	main:

		PREPROCESSINTERVALS (
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			params.bed,
			runname )


		COLLECTREADCOUNTS (
			bam,
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			PREPROCESSINTERVALS.out.interval_list )


		ANNOTATEINTERVALS (
			params.reference_fasta,
			params.reference_index,
			params.reference_dict,
			params.reference_gzi,
			PREPROCESSINTERVALS.out.interval_list )

		
		FILTERINTERVALS (
			PREPROCESSINTERVALS.out.interval_list,
			ANNOTATEINTERVALS.out.annotated_intervals,
			COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList() )

		
		INTERVALLISTTOOLS(
			FILTERINTERVALS.out.filter_intervals )


		DETERMINEGERMLINECONTIGPLOIDY(
			FILTERINTERVALS.out.filter_intervals,
			COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList(),
			params.contig_ploidy_priors )


		GERMLINECNVCALLER(
			DETERMINEGERMLINECONTIGPLOIDY.out.ploidy_calls,
			ANNOTATEINTERVALS.out.annotated_intervals,
			COLLECTREADCOUNTS.out.counts.collect().flatten().filter( ~/.*.counts.tsv$/ ).toList(),
			INTERVALLISTTOOLS.out.scatterout.flatten() )


		POSTPROCESSGERMLINECNVCALLS(
			DETERMINEGERMLINECONTIGPLOIDY.out.ploidy_calls,
			DETERMINEGERMLINECONTIGPLOIDY.out.sample_index_list.splitCsv(),
			GERMLINECNVCALLER.out.gcnv_dir.collect().flatten().filter( ~/.*calls$/ ).toList(),
			GERMLINECNVCALLER.out.gcnv_dir.collect().flatten().filter( ~/.*model$/ ).toList(),
			params.reference_dict )



	// emit:

}














// Aquí empezamos a encadenar los workflows: 



workflow {
	

	// Check that all parameters are OK
	CHECK_PARAMS ()
	

	// Download and concatenate samples
	if ( params.analysis.toUpperCase() =~ /[DM]/ ){
		
		DOWNLOAD(
			CHECK_PARAMS.out.controlsamples,
			CHECK_PARAMS.out.datasets )
	
	}


	// Mapping
	//YBQ: añadimos la opción de paralelizar el mapping 
	if ( params.analysis.toUpperCase().contains("M") ) {

		if ( params.parallel_mapping == true ){

			FASTP (DOWNLOAD.out.fastq) //process

			fastq_split=FASTP.out.reads //adaptamos el canal 
                //.map { sample, file -> [ file.getName().split('.')[0], file ]
                .map{ sample, file -> file }
                .collect()
                .flatten()
                //.fromFilePairs()
                //.map { file -> [ file.getParent(), file ]}
                //.map { file -> [ fromFilePairs("${file.getParent()}/*R{1,2}.fastp.fastq.gz") ]}
                .map { file -> [ file.getName().split('_')[0], file ] }
                .groupTuple()
                //.map{sample, file -> [sample.split('\\.')[1], file]} // con esto le puedo quitar el 0001 y solo dejarle el sample name.
                .flatten()
                .collate(3)
                //.view(),


			MAPPING( 
			fastq_split )
		
			//bam = MAPPING.out.bam

			bamstomerge= MAPPING.out.bam
                .map{sample, bam, bai -> [sample.split('\\.')[1], bam]} // con esto le puedo quitar el 0001 y solo dejarle el sample name.
                .groupTuple()
                //.view()

                bamstomerge.view()


                MERGEBAM(
                        bamstomerge,
                        params.assembly
                        )

                MERGEBAM.out.bam.view()

				
			bam = MERGEBAM.out.bam


		} else {

			MAPPING( 
			DOWNLOAD.out.fastq )
		
		bam = MAPPING.out.bam
		}

	}  


	if ( params.analysis.toUpperCase() =~ /[QSGCX]/ ){
		
		if ( params.analysis.toUpperCase().contains("M") ) {
					
			if ( params.parallel_mapping == true ){

				bam =  MERGEBAM.out.bam // si está en paralelo, asignamos a BAM la salida de MERGEBAM

			} else {

				bam = MAPPING.out.bam

			}
	
		} else {
			 // si no hay "M", comprobamos que haya un "BAM" local en la carpeta input. 
			LOCALBAM (
				params.input,
				CHECK_PARAMS.out.samples2analyce )
	
			bam = LOCALBAM.out.bam
		}
	}





	// Quality check
	if ( params.analysis.toUpperCase().contains("Q") ) {

		QUALITYCHECK(
			bam,
			CHECK_PARAMS.out.runname)
	
	}


/// mosdepth_bed -> si queremos el mosdepth bed que se genera para crear la base de datos
if ( params.mosdepth_bed == true ) {

	MOSDEPTH_COV(
			bam,
			params.bed,
			params.padding )
	
}

///keep el cram en una carpeta nueva llamada /cram
if ( params.keep_cram == true ) {

	BAM2CRAM (
		bam,
		params.reference_fasta,
		params.reference_index,
		params.reference_dict,
		params.reference_gzi,
		params.scratch ) //proceso
	
}


	// SNV calling
	if ( params.analysis.toUpperCase().contains("S") ) { 

		// Sample selection using JOIN function
		if (params.vc_tools.toLowerCase().split(',').contains('all')){

			SNVCALLING ( bam.join(CHECK_PARAMS.out.samples2analyce) )
			// ahora mismo, no se puede hacer en paralelo todas las herramientas, así que si 
			//vc_tools="all", ejecutamos el workflow normal (SNVCALLING)
			// A añadir: si split=yes, hacer en paralelo gatk y dragen y normal deepvariant y merge todo. 

		} else {

			if ( params.parallel_calling == true ) {
				SPLIT_BAM( 
				bam.join(CHECK_PARAMS.out.samples2analyce) 
				) //proceso
				bam = SPLIT_BAM.out.transpose()
				if ( params.vc_tools.toLowerCase().split(',').contains('gatk') ) {
					vcf = PARALLEL_GATKCALLING( bam ) //workflow
					vcf_gatk=PROCESS_GATKCALLING (vcf) //workflow
				} 
				
				if ( params.vc_tools.toLowerCase().split(',').contains('dragen') ){ 
					vcf = PARALLEL_DRAGENCALLING( bam ) //workflow
					vcf_dragen=PROCESS_DRAGENCALLING (vcf) //workflow
				} // YBQ: ¿FALTA UNIR AMBOS VCF FINALES EN EL CASO DE QUE SE EJECUTEN LOS DOS?
				//

			} else {
					bam = bam.join(CHECK_PARAMS.out.samples2analyce)
					// si no hacemos el parallel vcf calling: 

					if ( params.vc_tools.toLowerCase().split(',').contains('gatk') ) {
						HAPLOTYPECALLER (
								bam,
								params.bed,
								params.intervals,
								params.padding,
								params.reference_fasta,
								params.reference_index,
								params.reference_dict,
								params.reference_gzi,
								params.scratch ) //proceso
						vcf = HAPLOTYPECALLER.out.vcf
						vcf_gatk=PROCESS_GATKCALLING (vcf) //workflow
					}
					if ( params.vc_tools.toLowerCase().split(',').contains('dragen') ){ 
						STR_MODEL_DRAGEN(
							bam,
							params.reference_fasta,
							params.reference_index,
							params.reference_dict,
							params.reference_gzi,
							params.reference_str,
							params.scratch ) //proceso

						HAPLOTYPECALLER_DRAGEN(
							bam.join(STR_MODEL_DRAGEN.out.strmodel),
							params.bed,
							params.intervals,
							params.padding,
							params.reference_fasta,
							params.reference_index,
							params.reference_dict,
							params.reference_gzi,
							params.scratch ) //proceso
						vcf = HAPLOTYPECALLER_DRAGEN.out.vcf
						vcf_dragen=PROCESS_DRAGENCALLING (vcf) //workflow
					}
					if ( params.vc_tools.toLowerCase().split(',').contains('deepvariant') ) {
						DEEPVARIANT (
							bam.join(CHECK_PARAMS.out.samples2analyce),
							params.bed,
							params.intervals,
							params.reference_fasta,
							params.reference_index,
							params.reference_dict,
							params.reference_gzi,
							params.capture,
							params.scratch )

						vcf_deepvariant=FILTER_VCF_DEEPVARIANT (
							DEEPVARIANT.out.vcf,
							params.assembly,
							"deepvariant" )

					}
 


			}
		}

		
	
	}





	// SNV calling GVCF mode
	if ( params.analysis.toUpperCase().contains("G") ) {

		COMBINEDSNVCALLING ( 
			bam,
			CHECK_PARAMS.out.runname )
	}





	//////////////////// GUR: en SNV annotation tengo que editar esto para incluir K,E,T (añado else if, si no lo que ocurre es que chace dos checks_ LOCALBAM para el calling y LOCAL VCF para la anotacion y ambos checks los hace a la vez cuando todavia no ha corrido el calling entonces no hay bam, entonces con esto evitamos que compruebe el LOCALVCF cuando corre calling y anotacion)
	// SNV annotation
	if ( params.analysis.toUpperCase().contains("A") ) {
		if ( params.analysis.toUpperCase().contains("S") ){
			if ( params.vc_tools.toLowerCase().split(',').contains('all') ) {
			
			ANNOTATION ( SNVCALLING.out.finalvcf )

		} else if (params.vc_tools.toLowerCase().split(',').contains('gatk')) {
			
			ANNOTATION ( vcf_gatk )

		} else if ( params.vc_tools.toLowerCase().split(',').contains('dragen') ) {
			
			ANNOTATION ( vcf_dragen )

		} else if ( params.vc_tools.toLowerCase().split(',').contains('deepvariant') ) {
			
			ANNOTATION ( vcf_deepvariant )

		}
		} else { //si no viene del variant calling
			
			LOCALVCF (
				params.input,
				params.reference_fasta,
				CHECK_PARAMS.out.samples2analyce )

			vcf = LOCALVCF.out.vcf
			ANNOTATION ( vcf )
		}


		//ANNOTATION ( vcf )
	}





// yoli:EN LAS CNVS AÚN NO HEMOS AUTOMATIZADO EL PARALLEL MAPPING. 
	// CNV calling
	if ( params.analysis.toUpperCase().contains("C") ) {
		if ( params.analysis.toUpperCase().contains("M") ) {
			
			bam = MAPPING.out.bam			
			bam_collecteted = MAPPING.out.bam.collect().flatten().filter( ~/.*bam$/ ).toList()
			bai_collecteted = MAPPING.out.bam.collect().flatten().filter( ~/.*bai$/ ).toList()

		} else {

			LOCALBAM_CNV (
				params.input,
				CHECK_PARAMS.out.controlsamples )
			
			bam = LOCALBAM_CNV.out.bam
			bam_collecteted = LOCALBAM_CNV.out.bam.collect().flatten().filter( ~/.*bam$/ ).toList()
			bai_collecteted = LOCALBAM_CNV.out.bam.collect().flatten().filter( ~/.*bai$/ ).toList()

		}


		if ( params.capture.toUpperCase() == "P" ) {

			// CNV calling and annotation
			CNVCALLING ( 
				bam_collecteted,
				bai_collecteted, 
				CHECK_PARAMS.out.runname,
				CHECK_PARAMS.out.samples2analyce_file,
				bam )
		}


		if ( params.capture.toUpperCase() == "E" ) {

			// CNV calling and annotation
			CNVCALLING_WES ( 
				bam, 
				CHECK_PARAMS.out.runname )
		}


		if ( params.capture.toUpperCase() == "G" ) {

			// CNV calling and annotation
			CNVCALLING_WGS ( 
				bam_collecteted,
				bai_collecteted, 
				CHECK_PARAMS.out.runname )
		}


	}





}


















// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/13-1509.final.vcf --sample 13-1509 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/18-0439_filtered.vcf --sample 18-0439 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/18-2109_filtered.vcf --sample 18-2109 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results
// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --final_vcf /mnt/genetica4/gonzalo/muestras_ionut/21-0583.final.vcf --sample 21-0583 --analysis A --output /mnt/genetica4/gonzalo/muestras_ionut/results

// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --analysis M

// nextflow run ../annotation_nextflow/main.nf -profile hg19 --final_vcf head_21-4834-CES178.final.vcf --sample 21-4834-CES178 -with-report report.html -with-trace -with-timeline timeline.html -with-dag flowchart.png


// nextflow run /home/gonzalo/nextflowtest/annotation_nextflow/main.nf -profile hg19 --analysis A --input /mnt/genetica6/reanotacion/vcfs/ --output /mnt/genetica6/reanotacion/results/ -with-report report.7.html