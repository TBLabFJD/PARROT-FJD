




process BS_CHECK {
	label "basespace"
	publishDir "${params.output}/", mode: 'copy'

	input:
		val project
		val baseuser
		path samples
		val analysis

	output:
		path "projects.txt", emit: bsproyects
		path "controlsamples.txt", emit: controlsamples
		path "samples2analyce.txt", emit: samples2analyce
		path "datasets.txt", emit: datasets

	script:
		def baseuser_config = baseuser ? "--config ${baseuser} " : ''
		
		if(samples && analysis.contains("C"))
			"""
			# controlsamples: all samples need to be mapped for CNV calling 
			# samples2analyce: results (SNVs and CNVs) are only reported for the specified sample(s)

			# List all projects in the BaseSpace account
			bs list projects -f csv -F Name ${baseuser_config} > projects.txt

			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_config} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# Check that the specified samples exist inside the project 
			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi
			done


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""


		else if(samples) 
			"""
			# controlsamples: only the specifies sammple(s) is(are) mapped 
			# samples2analyce: results (SNVs) are only reported for the specified sample(s)

			# List all projects in the BaseSpace account
			bs list projects -f csv -F Name ${baseuser_config} > projects.txt


			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_config} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# Check that the specified samples exist inside the project
			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi   
			done

			# controlsamples are the same ones as samples2analyce
			cat samples2analyce.txt > controlsamples.txt


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""
		

		else 
			"""
			# controlsamples: all samples need to be mapped for SNV and CNV calling 
			# samples2analyce: results (SNVs and CNVs) are reported for all samples

			# List all projects in the BaseSpace account
			bs list projects -f csv -F Name ${baseuser_config} > projects.txt

			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName \
				${baseuser_config} --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# samples2analyce are the same ones as controlsamples
			cat controlsamples.txt > samples2analyce.txt


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""
		
}

process BS_CHECK_GUR {
	
	publishDir "${params.output}/", mode: 'copy'

	input:
		val project
		val baseuser
		path samples
		val analysis

	output:
		path "projects.txt", emit: bsproyects
		path "controlsamples.txt", emit: controlsamples
		path "samples2analyce.txt", emit: samples2analyce
		path "datasets.txt", emit: datasets

	script:
		if(samples && analysis.contains("C"))
			"""
			# controlsamples: all samples need to be mapped for CNV calling 
			# samples2analyce: results (SNVs and CNVs) are only reported for the specified sample(s)

			# List all projects in the BaseSpace account
			${baseuser}bs list projects -f csv -F Name > projects.txt

			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# Check that the specified samples exist inside the project 
			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi
			done


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			# edit el 21/05/2025: cambia la manera de descargar y ya no hay que hacer grep _L y tal, si se quieren descargar algunos antiguos igual si hace falta 
			${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
			# ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""


		else if(samples) 
			"""
			# controlsamples: only the specifies sammple(s) is(are) mapped 
			# samples2analyce: results (SNVs) are only reported for the specified sample(s)

			# List all projects in the BaseSpace account
			${baseuser}bs list projects -f csv -F Name > projects.txt


			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# Check that the specified samples exist inside the project
			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the project '${project}'\n"
					exit 1
				fi   
			done

			# controlsamples are the same ones as samples2analyce
			cat samples2analyce.txt > controlsamples.txt


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			# edit el 21/05/2025: cambia la manera de descargar y ya no hay que hacer grep _L y tal, si se quieren descargar algunos antiguos igual si hace falta 
			${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
			# ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""
		

		else 
			"""
			# controlsamples: all samples need to be mapped for SNV and CNV calling 
			# samples2analyce: results (SNVs and CNVs) are reported for all samples

			# List all projects in the BaseSpace account
			${baseuser}bs list projects -f csv -F Name > projects.txt

			# Check if the given project exist and if so, get all the sample names
			if grep -Fq ${project} projects.txt; then
				${baseuser}bs list biosample --sort-by=BioSampleName -f csv -F BioSampleName --project-name=${project} | sort | uniq > controlsamples.txt    
			else
				>&2 echo "ERROR: Project '${project}' does not exist in basespace\n"
				exit 1
			fi

			# samples2analyce are the same ones as controlsamples
			cat controlsamples.txt > samples2analyce.txt


			# Get the datasets id. 
			# First the last appsession Id is retrieved. 
			# Appsessions are the analysis done in a project. 
			# We assume that these are basecalling and that the last one is the correct one.
			appsession_id=\$(${baseuser}bs list appsession -f csv -F Id --project-name "${project}" | tail -n 1)

			# List all the Output.Datasets (folders containing the reads per sample and per lane) and 
			# filter to keep the ones containing the pattern *_L* to avoid duplicates.
			# edit el 21/05/2025: cambia la manera de descargar y ya no hay que hacer grep _L y tal, si se quieren descargar algunos antiguos igual si hace falta 
			${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name > datasets.txt
			# ${baseuser}bs appsession property get -i "\${appsession_id}" --property-name="Output.Datasets" -f csv -F Id -F Name | grep "_L" | grep -v "Undetermined" > datasets.txt
			"""
		
}











process LOCAL_CHECK {
	label "bioinfotools"
	publishDir "${params.output}/", mode: 'copy'

	input:
		path input
		path samples
		val analysis

	output:
		path "controlsamples.txt", emit: controlsamples
		path "samples2analyce.txt", emit: samples2analyce

	script:
		if ( analysis.contains("M") )      { extension_local_check = "(.fq|.fq.gz|.fastq|.fastq.gz)" }
		else if ( analysis.contains("Q") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("S") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("G") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("C") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("T") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("H") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("X") ) { extension_local_check = "(.bam|.cram)" }
		else if ( analysis.contains("A") ) { extension_local_check = "(.vcf|.vcf.gz)" }
		else if ( analysis.contains("N") ) { extension_local_check = "(.tsv|.bed)" }

		if(samples && analysis.contains("C"))
			"""
			ls ${input} | grep "${extension_local_check}\$" -P | sed 's/_.*//' | sed 's/\\..*//' | sed 's/${extension_local_check}\$//' -r | sort | uniq > controlsamples.txt
			

			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep -o \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the input folder '${input}'\n"
					exit 1
				fi
			done
			"""


		else if(samples) 
			"""
			ls ${input} | grep "${extension_local_check}\$" -P | sed 's/_.*//' | sed 's/\\..*//' | sed 's/${extension_local_check}\$//' -r | sort | uniq > controlsamples.txt
			

			for sample in \$(cat ${samples}); do
				if grep -q \${sample} controlsamples.txt; then
					grep -o \${sample} controlsamples.txt >> samples2analyce.txt
				else
					>&2 echo "ERROR: Sample '\${sample}' does not exist in the input folder '${input}'\n"
					exit 1
				fi
			done


			cat samples2analyce.txt > controlsamples.txt
			"""
		

		else 
			"""
			ls ${input} | grep "${extension_local_check}\$" -P | sed 's/_.*//' | sed 's/\\..*//' | sed 's/${extension_local_check}\$//' -r | sort | uniq > controlsamples.txt


			cat controlsamples.txt > samples2analyce.txt
			"""
		
}




// process BS_COPY {
// 	label "basespace"

// 	publishDir "${params.output}/fastq", mode: 'copy'
// 	maxRetries 4
// 	errorStrategy { task.attempt in 4 ? 'retry' : 'ignore' }


// 	input:
// 		val project
// 		val sample2download
// 		val baseuser

// 	output:
// 		tuple \
// 			val(sample2download_config), \
// 			path("${sample2download_config}_R1.fastq.gz"), \
// 			path("${sample2download_config}_R2.fastq.gz"), emit: fastq

// 	script:
// 		sample2download_config = sample2download[0]
// 		def baseuser_config = baseuser ? "--config ${baseuser} " : ''
// 		"""
// 		echo ${sample2download_config} > ${sample2download_config}.sample.txt
// 		bs download biosample ${baseuser_config} -q -n "${sample2download_config}" --exclude '*' --include '*.fastq.gz'


// 		cat ${sample2download_config}*/*_R1*fastq.gz > ${sample2download_config}_R1.fastq.gz
// 		cat ${sample2download_config}*/*_R2*fastq.gz > ${sample2download_config}_R2.fastq.gz
// 		"""
// }





process BS_COPY {
	label "basespace"

	publishDir "${params.output}/fastq", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }


	input:
		val project
		val sample2download
		val baseuser
		path datasets

	output:
		tuple \
			val(sample2download_config), \
			path("${sample2download_config}_R1.fastq.gz"), \
			path("${sample2download_config}_R2.fastq.gz"), emit: fastq

	script:
		sample2download_config = sample2download[0]
		def baseuser_config = baseuser ? "--config ${baseuser} " : ''
		"""
		grep "${sample2download_config}" ${datasets} | cut -d "," -f 1 | while read id; do
		echo "\${id}"
		bs download dataset -i "\${id}" --extension="fastq.gz"
		done

		cat *_R1*fastq.gz > ${sample2download_config}_R1.fastq.gz
		cat *_R2*fastq.gz > ${sample2download_config}_R2.fastq.gz
		"""
}

process BS_COPY_GUR {

	publishDir "${params.output}/fastq", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }


	input:
		val project
		val sample2download
		val baseuser
		path datasets

	output:
		tuple \
			val(sample2download_config), \
			path("${sample2download_config}_R1.fastq.gz"), \
			path("${sample2download_config}_R2.fastq.gz"), emit: fastq

	script:
		sample2download_config = sample2download[0]
		"""
		grep "${sample2download_config}" ${datasets} | cut -d "," -f 1 | while read id; do
		echo "\${id}"
		${baseuser}bs download dataset -i "\${id}" --extension="fastq.gz"
		done

		cat *_R1*fastq.gz > ${sample2download_config}_R1.fastq.gz
		cat *_R2*fastq.gz > ${sample2download_config}_R2.fastq.gz
		"""
}










process FASTQ_CONCATENATION {
	label "bioinfotools"

	input:
		path inputdir
		val sample2concat

	output:
		tuple \
			val(sample2concat_config), \
			path("${sample2concat_config}_R1.fastq.gz"), \
			path("${sample2concat_config}_R2.fastq.gz"), emit: fastq

	script:
		sample2concat_config = sample2concat[0]
		"""
		echo ${sample2concat_config}
		nR1="\$(ls ${inputdir}/${sample2concat_config}*R1*.f*q* | wc -l)"
		nR2="\$(ls ${inputdir}/${sample2concat_config}*R2*.f*q* | wc -l)"

		if [[ \$nR1 != \$nR2 ]]; then
			>&2 echo "Error: Number of R1 (\${nR1}) files do not match with the number of R2 (\${nR2}) files for the sample ${sample2concat_config}.\n"
			exit 1
		fi

		if [[ \$nR1 == 1 ]]; then
			ln -s ${inputdir}/${sample2concat_config}*R1*.f*q* ${sample2concat_config}_R1.fastq.gz
			ln -s ${inputdir}/${sample2concat_config}*R2*.f*q* ${sample2concat_config}_R2.fastq.gz
		elif [[ \$nR1 -gt 1 ]]; then
			cat ${inputdir}/${sample2concat_config}*R1*.f*q* > ${sample2concat_config}_R1.fastq.gz
			cat ${inputdir}/${sample2concat_config}*R2*.f*q* > ${sample2concat_config}_R2.fastq.gz 
		else
			>&2 echo "Error: wrong number of fastq. fastq files must match the pattern '*R[12]*.f*q*'\n"
			exit 1
		fi
		"""
}

		

process FASTP {
    label 'fastp'
    label "highcpu"
    label "highmem"

    input:
    tuple val(sample), path(forward), path(reverse)
        //val nsplit

    output:
                tuple val("${sample}"), \
                         path('000*.*_R*.fastp.fastq.gz'), emit: reads

    tuple val(sample), path('*.json')           , emit: json
    tuple val(sample), path('*.html')           , emit: html
    tuple val(sample), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(sample), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    script:

                """
                fastp \\
                --in1 ${forward} \\
                --in2 ${reverse} \\
                --out1 ${sample}_R1.fastp.fastq.gz \\
                --out2 ${sample}_R2.fastp.fastq.gz \\
                --json ${sample}.fastp.json \\
                --html ${sample}.fastp.html \\
                --thread 16 \\
                --detect_adapter_for_pe \\
                --disable_adapter_trimming \\
                --split 3

                """
}




process BWA {
	label "bwa"
	label "highcpu"
	label "highmem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(forward), path(reverse)
		path ref
		path bwa_amb
		path bwa_ann
		path bwa_pac
		path bwa_bwt
		path bwa_sa

	output:
		tuple \
			val(sample), \
			path("${sample}.mapped.bam"), emit: mapped_bam

	script:
		// sample  = fastq[0]
		// forward = fastq[1]
		// reverse = fastq[2]
		"""
		bwa mem -v 3 -t \$(nproc) -Y \\
		${ref} \\
		${forward} \\
		${reverse} |  samtools view -1 > ${sample}.mapped.bam

		"""
}





process FASTQTOSAM {
	label "gatk"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(forward), path(reverse)
		path scratch

	output:
		tuple \
			val(sample), \
			path("${sample}.unmapped.bam"), emit: unmapped_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_FastqToSam" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_FastqToSam" : ""

		"""

		${scratch_mkdir}

		gatk FastqToSam ${scratch_field} \
		--FASTQ ${forward} \
		--FASTQ2 ${reverse} \
		--OUTPUT ${sample}.unmapped.bam \
		--READ_GROUP_NAME ${sample} \
		--SAMPLE_NAME ${sample} \
		--PLATFORM illumina 
		"""
}




process MERGEBAMALIGNMENT {
	label "gatk"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(mapped_bam), path(unmapped_bam)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch

	output:
		tuple \
			val(sample), \
			path("${sample}.mapped.merged.bam"), emit: merged_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_MergeBamAlignment" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_MergeBamAlignment" : ""

		"""

		${scratch_mkdir}

		gatk MergeBamAlignment ${scratch_field} \
		--VALIDATION_STRINGENCY SILENT \
		--EXPECTED_ORIENTATIONS FR \
		--ATTRIBUTES_TO_RETAIN X0 \
		--ALIGNED_BAM ${mapped_bam} \
		--UNMAPPED_BAM ${unmapped_bam}  \
		--OUTPUT ${sample}.mapped.merged.bam \
		--REFERENCE_SEQUENCE ${ref} \
		--PAIRED_RUN true \
		--SORT_ORDER "unsorted" \
		--IS_BISULFITE_SEQUENCE false \
		--ALIGNED_READS_ONLY false \
		--CLIP_ADAPTERS false \
		--ADD_MATE_CIGAR true \
		--MAX_INSERTIONS_OR_DELETIONS -1 \
		--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
		--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
		--ALIGNER_PROPER_PAIR_FLAGS true \
		--UNMAP_CONTAMINANT_READS true
		"""
}









process MARKDUPLICATESSPARK {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	publishDir "${params.output}/mapping_stats", mode: 'copy', pattern: "marked_dup_metrics*"
	

	input:
		tuple val(sample), path(merged_bam)
		path scratch

	output:
		tuple \
			val(sample), \
			path("${sample}.dedupped.sorted.bam"), \
			path("${sample}.dedupped.sorted.bam.bai"), emit: deduppedsorted_bam
		tuple \
			val(sample), \
			path("marked_dup_metrics_${sample}.txt"), emit: dedupped_txt

	script:

		def scratch_field = scratch ? "--conf 'spark.local.dir=${scratch}/${sample}_MarkDuplicatesSpark'" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_MarkDuplicatesSpark" : ""

		"""
		${scratch_mkdir}

		gatk MarkDuplicatesSpark \
		-I ${merged_bam} \
		-O ${sample}.dedupped.sorted.bam \
		-M marked_dup_metrics_${sample}.txt \
		--remove-all-duplicates false \
		--optical-duplicate-pixel-distance 2500 \
		--read-validation-stringency SILENT \
		--create-output-bam-index true \
		${scratch_field}

		#chmod 777 \$(find . -user root) 
		chmod 777 *.dedupped.sorted.bam* .command.trace marked_dup_metrics*
		"""


}








process SORTSAM {
	label "gatk"	
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(merged_bam)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.sorted.bam"), \
			path("${sample}.sorted.bai"), emit: sorted_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_SortSam" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_SortSam" : ""
		"""

		${scratch_mkdir}

		gatk SortSam ${scratch_field} \
		--INPUT ${merged_bam} \
		--OUTPUT ${sample}.sorted.bam \
		--SORT_ORDER "coordinate" \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE false 
		"""
}








process SETTAGS {	
	label "gatk"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(sorted_bam), path(sorted_bai)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.tag.bam"), \
			path("${sample}.tag.bai"), emit: tagged_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_Settags" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_Settags" : ""
		
		"""

		${scratch_mkdir}

		gatk  SetNmMdAndUqTags ${scratch_field} \
		--INPUT ${sorted_bam} \
		--OUTPUT ${sample}.tag.bam \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE false \
		--REFERENCE_SEQUENCE ${ref}
		"""
}









process BASERECALIBRATOR {	
	label "gatk"
	publishDir "${params.output}/mapping_stats", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(deduppedsorted), path(deduppedsorted_bai)
		path ref
		path index
		path dict
		path reference_gzi
		path g1000_knownsites
		path g1000_knownsites_idx
		path mills_knownsites
		path mills_knownsites_idx
		path dbsnp_knownsites
		path dbsnp_knownsites_idx
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("before_recalibrated_bqsr_data_${sample}.recal.table"), emit: bqsr_table

	script:
		// def scratch_field = scratch ? "--TMP-DIR ${scratch}/${sample}_SortSam" : ""
		def g1000_knownsites_field = g1000_knownsites ? "--known-sites ${g1000_knownsites}" : ""
		def mills_knownsites_field = mills_knownsites ? "--known-sites ${mills_knownsites}" : ""
		def dbsnp_knownsites_field = dbsnp_knownsites ? "--known-sites ${dbsnp_knownsites}" : ""
		"""
		gatk BaseRecalibrator \
		-R ${ref} \
		-I ${deduppedsorted} \
		--use-original-qualities \
		${g1000_knownsites_field} \
		${mills_knownsites_field} \
		${dbsnp_knownsites_field} \
		-O before_recalibrated_bqsr_data_${sample}.recal.table
		"""
}









process APPLYBQSR {	
	label "gatk"
	publishDir "${params.output}/bams", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(deduppedsorted), path(deduppedsorted_bai), path(bqsr_table)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		val assembly
		
	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.bam"), \
			path("${sample}.${assembly}.bai"), emit: bam
		tuple \
			val(sample), \
			path("${sample}.${assembly}.bam.md5"), emit: md5 

	script:
		// def scratch_field = scratch ? "--TMP-DIR ${scratch}/${sample}_SortSam" : ""
		"""
		gatk ApplyBQSR \
		-R ${ref} \
		-I ${deduppedsorted}  \
		--bqsr ${bqsr_table} \
		-O ${sample}.${assembly}.bam \
		--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
		--add-output-sam-program-record \
		--create-output-bam-md5 \
		--use-original-qualities \
		--create-output-bam-index
		"""
}


process MERGEBAM{
        label "bioinfotools"
	publishDir "${params.output}/bams", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

        input:
                tuple val(sample), path(bam)
                val assembly

        output:
                tuple \
                        val(sample), \
                        path("${sample}.${assembly}.bam"),
						path("${sample}.${assembly}.bai"), emit: bam
                //tuple \
                //      val(sample), \
                //      path("${sample}.${assembly}.bam.md5"), emit: md5

        script:
                // def scratch_field = scratch ? "--TMP-DIR ${scratch}/${sample}_SortSam" : ""
                """
                samtools merge ${sample}.${assembly}_tmp.bam ${bam}
				samtools view -H ${sample}.${assembly}_tmp.bam  | sed "s/SM:[^\t]*/SM:${sample}/g" | samtools reheader - ${sample}.${assembly}_tmp.bam > ${sample}.${assembly}.bam
				rm ${sample}.${assembly}_tmp.bam
				samtools index ${sample}.${assembly}.bam ${sample}.${assembly}.bai
                """

}





process LOCALBAM {	
	label "bioinfotools"

	input:
		path inputdir
		val sample2analyce
		
	output:
		tuple \
			val(sample2analyce_config), \
			path("${sample2analyce_config}.bam"), \
			path("${sample2analyce_config}.bai"), emit: bam

	script:
		sample2analyce_config = sample2analyce[0]
		"""
		ln -s ${inputdir}/${sample2analyce_config}*.bam ${sample2analyce_config}.bam
		ln -s ${inputdir}/${sample2analyce_config}*.bai ${sample2analyce_config}.bai
		"""
}

// con ln -s crea un link simbólico

process LOCALCRAM {	
	label "bioinfotools"

	input:
		path inputdir
		val sample2analyce
		
	output:
		tuple \
			val(sample2analyce_config), \
			path("${sample2analyce_config}.cram"), \
			path("${sample2analyce_config}.crai"), emit: cram

	script:
		sample2analyce_config = sample2analyce[0]
		"""
		ln -s ${inputdir}/${sample2analyce_config}*.cram ${sample2analyce_config}.cram
		ln -s ${inputdir}/${sample2analyce_config}*.crai ${sample2analyce_config}.crai
		"""
}










process MOSDEPTH {
	label "bioinfotools"
	publishDir "${params.output}/qc/mosdepth/", mode: 'copy', pattern: "*txt"
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed

	output:
		path("${sample}_qc.*"), emit: mosdepth

	script:
		def bed_field = bed ? "--by ${bed}" : ""	
		
		"""
		mosdepth -x --no-per-base ${bed_field} ${sample}_qc ${bam}
		"""
}


/*
process MOSDEPTH_JOIN {
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path ""
		val readthreshold
		val analysis
		val run

	output:
		path("minCovFilterResults_${analysis}.txt"), emit: coverage
		path("discarded_samples_${analysis}.txt"), emit: discarded_samples optional true 


	script:
		
		"""
		for file in *_qc.mosdepth.region.dist.txt; do

			sample="\$(basename \${file} _qc.mosdepth.region.dist.txt)"

			cov=\$(awk -v reads="${readthreshold}" '{if(\$1=="total" && \$2==reads){print \$3}}' \${file})
			pass=\$(awk 'BEGIN{if ('\${cov}'>'0.9') print 0}')

			if [ "\${pass}" = "0" ]; then
				echo -e \${sample}"\\t"${readthreshold}"\\t"\${cov}"\\tPASSED\\t"${run} >> minCovFilterResults_${analysis}.txt
			else
				echo -e \${sample}"\\t"${readthreshold}"\\t"\${cov}"\\tFAILED\\t"${run} >> minCovFilterResults_${analysis}.txt
				echo -e \${sample} >> discarded_samples_${analysis}.txt
			fi
		done
		"""
}

*/

process MOSDEPTH_PLOT {
	label "bioinfotools"
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path ""
		path projectDir

	output:
		path("mosdepth.region.dist.html"), emit: plot
		

	script:
		
		"""
		python ${projectDir}/tasks/plot-dist.py *.mosdepth.region.dist.txt -o mosdepth.region.dist.html
		"""
}





process MOSDEPTH_COV {
	label "bioinfotools"
	publishDir "${params.output}/qc/mosdepth_cov", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val padding

	output:
		tuple \
			val(sample), \
			path("${sample}*bed"), emit: mosdepth

	script:
		
/*		if(bed)
			"""
			mosdepth --quantize 10: -n -x ${sample}_cov ${bam}
			zcat ${sample}_cov.quantized.bed.gz > ${sample}.global.quantized.bed

			awk -v pad="${padding}" '{print \$1"\\t"\$2-pad"\\t"\$3+pad"\\t"\$4}' ${bed} | bedtools merge -i - > bedpadding.bed
			bedtools intersect -a bedpadding.bed -b ${sample}_cov.quantized.bed.gz -wb | awk '{print \$1"\\t"\$2"\\t"\$3"\\t"\$NF}' > ${sample}.padding.quantized.bed

			
			#bedtools coverage -a BED -b BAM -mean -sorted > MeanCoverageBED.bedgraph
			#bedtools map -a BED -b BAM -mean -sorted > MeanCoverageBED.bedgraph
			"""

		else*/
			"""
			mosdepth --quantize 10: -n -x ${sample}_cov ${bam}
			zcat ${sample}_cov.quantized.bed.gz > ${sample}.${params.technique}.v41.${params.current_date}.global.quantized.bed
			"""
}








process GENOMECOV {
	label "bioinfotools"
	publishDir "${params.output}/qc/genomecov/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed

	output:
		tuple \
			val(sample), \
			path("${sample}*bed"), emit: genomecov

	script:
		
		if(bed)
			"""
			samtools view -b -F 2304 ${bam} | bedtools genomecov -bga -ibam stdin > ${sample}.genomecov.bed
			#bedtools intersect -a ${sample}.genomecov.bed -b ${bed} > ${sample}.panelcov.bed
			bedtools intersect -a ${sample}.genomecov.bed -b ${bed} -wb > ${sample}.panelcov.bed
			"""

		else
			"""
			samtools view -b -F 2304 ${bam} | bedtools genomecov -bga -ibam stdin > ${sample}.genomecov.bed
			"""
}
// 3332
// read unmapped (0x4)4
// not primary alignment (0x100) 256
// read is PCR or optical duplicate (0x400) 1024
// supplementary alignment (0x800) 2048

// 2304
// not primary alignment (0x100) 256
// supplementary alignment (0x800) 2048


process SAMTOOLS_FLAGSTAT {
	label "bioinfotools"
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.samtools_flagstat.txt"), emit: flagstat

	script:

		"""
		samtools flagstat ${bam} > ${sample}.samtools_flagstat.txt
		"""
}






process READ_LENGTH_STATS {
	label "bioinfotools"
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.read_lenghts.txt"), emit: read_length_stats

	script:

		"""
		samtools view -F 2304 ${bam} | cut -f 10 | perl -ne 'chomp;print length(\$_) . "\\n"' > ${sample}.read_lenghts.txt
		"""
}



process SEQUENCING_QUALITY_SCORES {
	label "bioinfotools"
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.quality.txt"), emit: quality

	script:

		"""
		samtools view -F 2304 ${bam} | cut -f 11 | grep -o . | awk '{c[\$0]++}END{for(l in c){print c[l], l}}' | sort -n > ${sample}.quality.txt
		"""
}



process SEQUENCING_CG_AT_CONTENT {
	label "bioinfotools"
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.CG_AT.txt"), emit: cg_at

	script:

		"""
		samtools view -F 2304 ${bam} | cut -f 10 | grep -o . | awk '{c[\$0]++}END{for(l in c){print c[l], l}}' | sort -n > ${sample}.CG_AT.txt
		"""
}




process NREADS_NONDUP_UNIQ {
	label "bioinfotools"
	// publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.nreads_nondup_uniq.txt"), emit: nreads_nondup_uniq

	script:

		"""
		samtools view -F 3332 ${bam} | wc -l > ${sample}.nreads_nondup_uniq.txt
		"""
}



process QC_SUMMARY {
	label "bioinfotools"
	label "highmem"
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(""), path("${sample}.samtools_flagstat.txt"), path("${sample}.read_lenghts.txt"), path("${sample}.quality.txt"), path("${sample}.CG_AT.txt"), path("${sample}.nreads_nondup_uniq.txt")
		path bed
		path projectDir

	output:
		tuple \
			val(sample), \
			path("${sample}*txt"), emit: quality_summary

	script:
		def bed_path = bed ? "--bed_path ${bed} " : ''
		def panelcov_path = bed ? "--panelcov_path ${sample}.panelcov.bed " : ''
		def bedoutpath = bed ? "--bedout ${sample}.library.stats.txt " : ''

		"""
		nreads_nondup_uniq="\$(cat ${sample}.nreads_nondup_uniq.txt)"

		Rscript ${projectDir}/tasks/quality_summary.R \\
		--samplename ${sample} \\
		--genomecov_path ${sample}.genomecov.bed \\
		--quality_path ${sample}.quality.txt \\
		--samtools_flagstat_path ${sample}.samtools_flagstat.txt \\
		--read_lenghts_path ${sample}.read_lenghts.txt \\
		--cg_at_path ${sample}.CG_AT.txt \\
		${panelcov_path}\\
		${bed_path}\\
		${bedoutpath}\\
		--output ${sample}.quality.summary.txt \\
		--nreads_nondup_uniq \${nreads_nondup_uniq}
		"""
}




process RUN_QC_CAT {
	label "bioinfotools"
	publishDir "${params.output}/qc/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path("")
		val runname

	output:
		tuple \
			val(runname), \
			path("${runname}.quality.summary.txt"), emit: quality_summary

	script:
		"""
		first_file="\$(ls *quality.summary.txt | head -n 1)"
		head -n 1 \${first_file} > ${runname}.quality.summary.txt.tmp
		tail -n 1 -q *quality.summary.txt >> ${runname}.quality.summary.txt.tmp
		mv ${runname}.quality.summary.txt.tmp ${runname}.quality.summary.txt
		"""
}



//params.chroms = 'chr{X,Y}'
params.chroms = 'chr{{1..22},X,Y}'

process SPLIT_BAM {
	label "bioinfotools"
	label "highcpu"
	label "highmem"
	tag { sample }
	//publishDir "${params.output}/split_bams", mode: 'copy'
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    input:
    tuple val(sample), path(bam), path(bai)

    output:
		tuple \
			val(sample), \
			path("parallel_bams/*.bam"), emit: parallel_bams

    //path("parallel_bams/*.bam"), emit: parallel_bams

    script:
    """
	mkdir parallel_bams
    for chrom in ${params.chroms}
    do
        samtools view \\
            -o "parallel_bams/${bam.baseName}.\${chrom}.bam" \\
            "${bam}" \\
            "\${chrom}"
		samtools index "parallel_bams/${bam.baseName}.\${chrom}.bam"
    done
    """
}



process PARALLEL_HAPLOTYPECALLER {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	//publishDir "${params.output}/", mode: 'copy'
	//publishDir "${params.output}/out_parallel_vcfs/", mode: 'copy'
	tag { bam }

	input:
		tuple val(sample), path(bam, stageAs: 'parallel_bams/*')
		//tuple val(sample), path(bam, stageAs: 'parallel_bams/*')
		//path bam, stageAs: 'parallel_bams/*'
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		//path("${bam.baseName}.vcf"), emit: vcf
		//path("${bam.baseName}.vcf.idx"), emit: vcf_idx
		tuple \
			val(sample), \
			path("${bam.baseName}.vcf.idx"), emit: vcf_idx
		tuple \
			val(sample), \
			path("${bam.baseName}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${bam.baseName}_HaplotypeCaller" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${bam.baseName}_HaplotypeCaller" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}
		mkdir parallel_vcfs
		gatk --java-options "-Xmx${params.mediummem}g" \
		HaplotypeCaller ${scratch_field} \
		-R ${ref} \
		-I ${bam} \
		-O "${bam.baseName}.vcf" \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}




process HAPLOTYPECALLER {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.vcf"), \
			path("${sample}.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCaller" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCaller" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk --java-options "-Xmx${params.mediummem}g" \
		HaplotypeCaller ${scratch_field} \
		-R ${ref} \
		-I ${bam} \
		-O ${sample}.vcf \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}


process MERGE_SPLIT_VCF {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	tag { sample }
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	//publishDir "${params.output}/parallel_vcfs", mode: 'copy'
	input:
		//path my_vcfs
		tuple val(sample), path(my_vcfs)
		//, stageAs: "${params.output}/out_parallel_vcfs"
		//path my_vcfs, stageAs: 'parallel_vcfs/*'
		path ref
		path scratch
		val program
		
	output:
		
		//path("merged.vcf"), emit: vcf
		tuple \
			val(sample), \
			path("${sample}.${program}.vcf"), \
			path("${sample}.${program}.vcf.idx"), emit: vcf
		//path("${my_vcfs.SimpleName}.vcf.idx"), emit: vcf_idx


	script:
		//def scratch_field   = scratch ? "--TMP_DIR ${scratch}/${my_vcfs.SimpleName}_mergesplitvcf" : ""	
		//def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${my_vcfs.SimpleName}_mergesplitvcf" : ""

		"""
		echo ${my_vcfs} | tr ' ' '\n' > vcfs.list
		gatk MergeVcfs \
		-R ${ref} \
		-I vcfs.list \
		-O "${sample}.${program}.vcf"
		"""
}



process SELECT_SNV {
	label "gatk"	
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }

	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.snp.vcf"), \
			path("${sample}.snp.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectsnvs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectsnvs" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include SNP \
		-O ${sample}.snp.vcf
		"""
}





process SELECT_INDEL {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }

	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.indel.vcf"), \
			path("${sample}.indel.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectindel" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectindel" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include INDEL \
		-O ${sample}.indel.vcf
		"""
}




process SELECT_MIX {
	label "gatk"	
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }
	input:
		tuple val(sample), path(vcf), path(idx)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.mix.vcf"), \
			path("${sample}.mix.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_selectmix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_selectmix" : ""

		"""
		${scratch_mkdir}

		gatk SelectVariants ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--select-type-to-include MIXED \
		--select-type-to-include MNP \
		--select-type-to-include SYMBOLIC \
		-O ${sample}.mix.vcf
		"""
}





process FILTRATION_SNV {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.sorted.snv.vcf"), \
			path("${sample}.filtered.sorted.snv.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filtersnvs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtersnvs" : ""
		def scratch_field2   = scratch ? "--TMP_DIR ${scratch}/${sample}_filtersnvs" : ""	

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "SOR > 3.0" --filter-name "SOR3" \
		-filter "FS > 60.0" --filter-name "FS60" \
		-filter "MQ < 40.0" --filter-name "MQ40" \
		-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
		-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
		-O ${sample}.filtered.snv.vcf

		gatk SortVcf ${scratch_field2} \
		-I ${sample}.filtered.snv.vcf \
		-O ${sample}.filtered.sorted.snv.vcf

		"""
}







process FILTRATION_INDEL {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.sorted.indel.vcf"), \
			path("${sample}.filtered.sorted.indel.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filterindel" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filterindel" : ""
		def scratch_field2   = scratch ? "--TMP_DIR ${scratch}/${sample}_filterindel" : ""	

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-filter "FS > 200.0" --filter-name "FS200" \
		-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
		-O ${sample}.filtered.indel.vcf

		gatk SortVcf ${scratch_field2} \
		-I ${sample}.filtered.indel.vcf \
		-O ${sample}.filtered.sorted.indel.vcf

		"""
}








process FILTRATION_MIX {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	tag { sample }

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.filtered.sorted.mix.vcf"), \
			path("${sample}.filtered.sorted.mix.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filtermix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtermix" : ""
		def scratch_field2   = scratch ? "--TMP_DIR ${scratch}/${sample}_filtermix" : ""	

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		-filter "QD < 2.0" --filter-name "QD2" \
		-filter "QUAL < 30.0" --filter-name "QUAL30" \
		-O ${sample}.filtered.mix.vcf 
		
		gatk SortVcf ${scratch_field2} \
		-I ${sample}.filtered.mix.vcf \
		-O ${sample}.filtered.sorted.mix.vcf
		
		"""
}






process MERGE_VCF {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	//publishDir "${params.output}/snvs", mode: 'copy'
	tag { sample }

	input:
		tuple val(sample), path(vcf_snv), path(idx_snv), path(vcf_indel), path(idx_indel), path(vcf_mix), path(idx_mix)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.gatkLabeled.vcf"), \
			path("${sample}.gatkLabeled.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--TMP_DIR ${scratch}/${sample}_filtermix" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filtermix" : ""

		"""
		${scratch_mkdir}

		gatk MergeVcfs ${scratch_field} \
		-R ${ref} \
		-I ${vcf_snv} \
		-I ${vcf_indel} \
		-I ${vcf_mix} \
		-O ${sample}.gatkLabeled.vcf
		"""
}



// process FILTER_VCF {	
// 	label "bioinfotools"
// 	publishDir "${params.output}/snvs", mode: 'copy'
	
// 	input:
// 		tuple val(sample), path(vcf_snv), path(idx_snv)
		
// 	output:
// 		tuple \
// 			val(sample), \
// 			path("${sample}.final.gatk.vcf.gz"),
// 			path("${sample}.final.gatk.vcf.gz.tbi"), emit: vcf

// 	script:

// 		"""
// 		bcftools filter -O z -o ${sample}.final.gatk.vcf.gz -i FILTER~"PASS" \
// 		-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
// 		${vcf_snv}

// 		tabix -p vcf ${sample}.final.gatk.vcf.gz
// 		"""
// }











process DEEPVARIANT {
	label "deepvariant"	
	label "highcpu"	
	label "highmem"	
	errorStrategy 'retry'
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val intervals
		path ref
		path index
		path dict
		path reference_gzi
		val capture
		path scratch

		
	output:
		tuple \
			val(sample), \
			path("${sample}.deepvariantLabeled.vcf.gz"), \
			path("${sample}.deepvariantLabeled.vcf.gz.tbi"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.deepvariantLabeled.gvcf.gz"), \
			path("${sample}.deepvariantLabeled.gvcf.gz.tbi"), emit: gvcf

	script:
		def capture_field   = capture != "G" ? "WES" : "WGS"
		def intervals_field = intervals ? "--regions ${bed}" : ""
		def scratch_field   = scratch ? "--intermediate_results_dir ${scratch}/${sample}_deepvariant" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_deepvariant" : ""

		"""
		${scratch_mkdir}

		run_deepvariant ${scratch_field} \
		--model_type=${capture_field} \
		--ref=${ref} \
		--reads=${bam} \
		--output_vcf=${sample}.deepvariantLabeled.vcf.gz \
		--output_gvcf=${sample}.deepvariantLabeled.gvcf.gz \
		--num_shards=\$(nproc) \
		${intervals_field} \
		"""
}




// process FILTER_VCF_DEEPVARIANT {	
// 	label "bioinfotools"
// 	publishDir "${params.output}/snvs", mode: 'copy'
	
// 	input:
// 		tuple val(sample), path(vcf_snv), path(idx_snv)
		
// 	output:
// 		tuple \
// 			val(sample), \
// 			path("${sample}.final.deepvariant.vcf.gz"),
// 			path("${sample}.final.deepvariant.vcf.gz.tbi"), emit: vcf

// 	script:

// 		"""
// 		bcftools filter -O z -o ${sample}.final.deepvariant.vcf.gz -i FILTER~"PASS" \
// 		-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
// 		${vcf_snv}

// 		tabix -p vcf ${sample}.final.deepvariant.vcf.gz
// 		"""
// }


///// (graci) MI PROCESO: PARALLEL DRAGEN: (STR + haplotype caller) -> paraleliza y ademas junta los dos procesos de STR y haplotype caller
process PARALLEL_HAPLOTYPECALLER_DRAGEN {	
	label "gatk"
	errorStrategy 'ignore'
	// publishDir "${params.output}/", mode: 'copy'
	//INFO: bam.baseName (sample.chromosome) -> file name without its extension: 23-0136.chr5
	//INFO: bam.SimpleName (cuando no hay cromosoma: sample ) -> file name without any extension: 23-0136
	// 

	input:
		tuple val(sample), path(bam, stageAs: 'parallel_bams/*')
		path ref
		path index
		path dict
		path reference_gzi
		path reference_str
		path scratch
		path bed
		val intervals
		val padding
		
	output:
		tuple \
			val(sample), \
			path("${bam.baseName}.vcf.idx"), emit: vcf_idx
		tuple \
			val(sample), \
			path("${bam.baseName}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${bam.BaseName}_mixmodeldragen" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${bam.BaseName}_mixmodeldragen" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk CalibrateDragstrModel ${scratch_field} \
    	-R ${ref} \
    	-I ${bam} \
    	-str ${reference_str} \
    	-O ${bam.baseName}_dragstr_model.txt

		gatk --java-options "-Xmx${params.mediummem}g" \
		HaplotypeCaller ${scratch_field} \
		--dragen-mode \
		--dragstr-params-path ${bam.baseName}_dragstr_model.txt \
		-R ${ref} \
		-I ${bam} \
		-O ${bam.baseName}.vcf \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}

		"""
}





process STR_MODEL_DRAGEN {	
	label "gatk"
	errorStrategy 'ignore'
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai)
		path ref
		path index
		path dict
		path reference_gzi
		path reference_str
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}_dragstr_model.txt"), emit: strmodel

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_strmodeldragen" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_strmodeldragen" : ""

		"""
		${scratch_mkdir}

		gatk CalibrateDragstrModel ${scratch_field} \
    	-R ${ref} \
    	-I ${bam} \
    	-str ${reference_str} \
    	-O ${sample}_dragstr_model.txt

		"""
}



process HAPLOTYPECALLER_DRAGEN {	
	label "gatk"
	errorStrategy 'ignore'
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai), path(str)
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.vcf"), \
			path("${sample}.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCallerDragen" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCallerDragen" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk --java-options "-Xmx${params.mediummem}g" \
		HaplotypeCaller ${scratch_field} \
		--dragen-mode \
		--dragstr-params-path ${str} \
		-R ${ref} \
		-I ${bam} \
		-O ${sample}.vcf \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}



process FILTRATION_DRAGEN {	
	label "gatk"
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(vcf), path(idx)
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.dragenLabeled.vcf"), \
			path("${sample}.dragenLabeled.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_filterdragen" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_filterdragen" : ""

		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-V ${vcf} \
		--filter-expression "QUAL < 10.4139" \
      	--filter-name "DRAGENHardQUAL" \
		-O ${sample}.dragenLabeled.vcf
		"""
}



// process FILTER_VCF_DRAGEN {	
// 	label "bioinfotools"

// 	publishDir "${params.output}/snvs", mode: 'copy'
// 	input:
// 		tuple val(sample), path(vcf_snv), path(idx_snv)
		
// 	output:
// 		tuple \
// 			val(sample), \
// 			path("${sample}.final.dragen.vcf.gz"),
// 			path("${sample}.final.dragen.vcf.gz.tbi"), emit: vcf

// 	script:

// 		"""
// 		bcftools filter -O z -o ${sample}.final.dragen.vcf.gz -i FILTER~"PASS" \
// 		-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
// 		${vcf_snv}

// 		tabix -p vcf ${sample}.final.dragen.vcf.gz
// 		"""
// }


process FILTER_VCF {	
	label "bioinfotools"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	publishDir "${params.output}/individual_callers_snvs", mode: 'copy'
	input:
		tuple val(sample), path(vcf_snv), path(idx_snv)
		val assembly
		val program
		
	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.${program}.vcf.gz"),
			path("${sample}.${assembly}.${program}.vcf.gz.tbi"), emit: vcf

	script:

		"""
		bcftools view -O z -o tmp.vcf.gz ${vcf_snv} 
		tabix -p vcf tmp.vcf.gz

		bcftools filter -O z -o ${sample}.${assembly}.${program}.vcf.gz -i 'FILTER~"PASS"' \
		-r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
		tmp.vcf.gz

		tabix -p vcf ${sample}.${assembly}.${program}.vcf.gz
		"""
}

//// GUR: PROCESO NUEVO PARA OBTENER EL VCF FINAL EN LA CARPETA SNVS (EL FINAL VCF)
process FINAL_VCF {	
	label "bioinfotools"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(vcf_snv), path(idx_snv)
		//tuple path(vcf_snv), path(idx_snv)
		val assembly
		val program
		
	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.${program}.${params.technique}.v41.${params.current_date}.final.vcf.gz"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.${assembly}.${program}.${params.technique}.v41.${params.current_date}.final.vcf.gz.tbi"), emit: index

	script:

		"""
		#convertir el vcf individual en el final y crearle su index (basicamente renombrarlo)
	    	cp ${sample}.${assembly}.${program}.vcf.gz ${sample}.${assembly}.${program}.${params.technique}.v41.${params.current_date}.final.vcf.gz
		tabix -p vcf ${sample}.${assembly}.${program}.${params.technique}.v41.${params.current_date}.final.vcf.gz
	
		"""
}

/////////// PROCESO CONVERSION BAM/CRAM
process BAM2CRAM {	
	label "bioinfotools"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	publishDir "${params.output}/cram", mode: 'copy'
	tag { sample }

	input:
		tuple val(sample), path(bam), path(bai)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${bam.baseName}.${params.technique}.v41.${params.current_date}.cram.crai"), emit: cram_idx
		tuple \
			val(sample), \
			path("${bam.baseName}.${params.technique}.v41.${params.current_date}.cram"), emit: cram

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_bam2cram" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_bam2cram" : ""

		"""
		${scratch_mkdir}
		samtools view -C -T ${ref} -o ${bam.baseName}.${params.technique}.v41.${params.current_date}.cram ${bam} 
		samtools index ${bam.baseName}.${params.technique}.v41.${params.current_date}.cram ${bam.baseName}.${params.technique}.v41.${params.current_date}.cram.crai
		"""
}

process CRAM2BAM {	
	label "bioinfotools"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	// publishDir "${params.output}/", mode: 'copy'

	input:
		tuple val(sample), path(cram), path(crai)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:

		tuple \
			val(sample), \
			path("${cram.baseName}.bam"), \
			path("${cram.baseName}.bai"), emit: bam

		//tuple \
		//	val(sample), \
		//	path("${bam.baseName}.bam.bai"), emit: bam_idx
		//tuple \
		//	val(sample), \
		//	path("${bam.baseName}.bam"), emit: bam

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_cram2bam" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_cram2bam" : ""

		"""
		${scratch_mkdir}
		samtools view -b -T ${ref} -o ${cram.baseName}.bam ${cram} 
		samtools index ${cram.baseName}.bam ${cram.baseName}.bai

		"""
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// GUR: PROCESOS NUEVOS FINALES PARA QUE PUEDA USARLO PARA ANOTAR -> se crea la carpeta snvs igual que en el MERGE_VCF_CALLERS Y el mismo .vcf de individual callers pasa a ser el archivo .final.vcf ///////

/*
process FINAL_GATK{	
	label "bioinfotools"
	errorStrategy 'ignore'

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(gatk_vcf), path(gatk_tbi)
		val assembly
		val program

	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz.tbi"), emit: index

	script:

		"""
		#convertir el vcf individual en el final y crearle su index (basicamente renombrarlo)
	    cp ${sample}.${assembly}.${program}.vcf.gz ${sample}.${assembly}.final.vcf.gz
		tabix -p vcf ${sample}.${assembly}.final.vcf.gz
		"""
}


process FINAL_DEEPVARIANT{	
	label "bioinfotools"
	errorStrategy 'ignore'

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(deepvariant_vcf), path(deepvariant_tbi)
		val assembly
		val program

	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz.tbi"), emit: index

	script:

		"""
		#convertir el vcf individual en el final y crearle su index (basicamente renombrarlo)
	    cp ${sample}.${assembly}.${program}.vcf.gz ${sample}.${assembly}.final.vcf.gz
		tabix -p vcf ${sample}.${assembly}.final.vcf.gz
		"""
}


process FINAL_DRAGEN{	
	label "bioinfotools"
	errorStrategy 'ignore'

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(dragen_vcf), path(dragen_tbi)
		val assembly
		val program

	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz.tbi"), emit: index

	script:

		"""
		#convertir el vcf individual en el final y crearle su index (basicamente renombrarlo)
	    cp ${sample}.${assembly}.${program}.vcf.gz ${sample}.${assembly}.final.vcf.gz
		tabix -p vcf ${sample}.${assembly}.final.vcf.gz
		"""
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



process MERGE_VCF_CALLERS {	
	label "bioinfotools"
	errorStrategy 'ignore'

	publishDir "${params.output}/snvs", mode: 'copy'
	input:
		tuple val(sample), path(gatk_vcf), path(gatk_tbi), path(deepvariant_vcf), path(deepvariant_tbi), path(dragen_vcf), path(dragen_tbi)
		val assembly
		path ref
		path projectDir

	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz"), emit: vcf

		tuple \
			val(sample), \
			path("${sample}.${assembly}.final.vcf.gz.tbi"), emit: index

	script:

		"""
		# Add program as sufix to sample names
		bcftools query -l ${deepvariant_vcf} > samples.deepvariant.txt
		sed -e 's/\$/.DV/' -i samples.deepvariant.txt
		bcftools reheader --samples samples.deepvariant.txt -o ${sample}.renamed.deepvariant.vcf.gz ${deepvariant_vcf}
		tabix -p vcf ${sample}.renamed.deepvariant.vcf.gz

		bcftools query -l ${dragen_vcf} > samples.dragen.txt
		sed -e 's/\$/.DR/' -i samples.dragen.txt
		bcftools reheader --samples samples.dragen.txt -o ${sample}.renamed.dragen.vcf.gz ${dragen_vcf}
		tabix -p vcf ${sample}.renamed.dragen.vcf.gz

		bcftools query -l ${gatk_vcf} > samples.gatk.txt
		sed -e 's/\$/.GK/' -i samples.gatk.txt
		bcftools reheader --samples samples.gatk.txt -o ${sample}.renamed.gatk.vcf.gz ${gatk_vcf}
		tabix -p vcf ${sample}.renamed.gatk.vcf.gz



		# Split multiallelic reacords as biallelic records and keep only records wth a minimum allele count of 1
		bcftools norm -m - -c s -f ${ref} ${sample}.renamed.deepvariant.vcf.gz | bcftools view --min-ac=1 -O z -o ${sample}.biallelic.deepvariant.vcf.gz 
		tabix -p vcf ${sample}.biallelic.deepvariant.vcf.gz
		bcftools norm -m - -c s -f ${ref} ${sample}.renamed.dragen.vcf.gz | bcftools view --min-ac=1 -O z -o ${sample}.biallelic.dragen.vcf.gz 
		tabix -p vcf ${sample}.biallelic.dragen.vcf.gz
		bcftools norm -m - -c s -f ${ref} ${sample}.renamed.gatk.vcf.gz | bcftools view --min-ac=1 -O z -o ${sample}.biallelic.gatk.vcf.gz 
		tabix -p vcf ${sample}.biallelic.gatk.vcf.gz



		# Merge VCFs from the 3 programs
		bcftools merge -O z -o ${sample}.merged.vcf.gz -m none ${sample}.biallelic.deepvariant.vcf.gz ${sample}.biallelic.dragen.vcf.gz ${sample}.biallelic.gatk.vcf.gz
		tabix -p vcf ${sample}.merged.vcf.gz



		# Get genotype, depth and variant allele depth
		bcftools query -f '[\\t%GT]\\n' ${sample}.merged.vcf.gz | sed 's/\\t//1' | sed 's/\\.\\/\\.//g' | sed 's/1\\/0/0\\/1/g' > GT.txt
		bcftools query -f '[\\t%DP]\\n' ${sample}.merged.vcf.gz | sed 's/\\t//1' | sed 's/\\.//g' > DP.txt
		bcftools query -f '[\\t%AD{1}]\\n' ${sample}.merged.vcf.gz | sed 's/\\t//1' | sed 's/\\.//g' > VD.txt
		bcftools query -f '[\\t%AD{0}]\\n' ${sample}.merged.vcf.gz | sed 's/\\t//1' | sed 's/\\.//g' > RD.txt



		# Calculate the average depth and variant allele depth
		awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' DP.txt > DP_mean.txt
		awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' VD.txt > VD_mean.txt
		awk -v OFMT=%.0f '{sum = 0; for (i = 1; i <= NF; i++) sum += \$i; sum /= NF; print sum}' RD.txt > RD_mean.txt
		paste -d, RD_mean.txt VD_mean.txt > AD_mean.txt

		# Get consensus genotype
		python ${projectDir}/tasks/consensus_GT.py GT.txt GT_consensus.txt GT_discordances.txt


		# Calculate variant allele depth (VAD)
		paste VD_mean.txt DP_mean.txt | awk -v OFMT=%.2f '{print(\$1/\$2)}' > VAF.txt




		# Software information
		cut -f 1 GT.txt | sed 's/.+/DV/' -r > DV.txt
		cut -f 2 GT.txt | sed 's/.+/DG/' -r > DG.txt
		cut -f 3 GT.txt | sed 's/.+/GK/' -r > GK.txt

		paste -d "\\t" DV.txt DG.txt GK.txt | perl -alne 'print join "_", @F' > SF.txt



		# Create sample information
		bcftools query -f '[:%GT:%DP:%AD{1}]\\n' ${sample}.merged.vcf.gz > FORMAT_SF.txt
		paste -d ":" GT_consensus.txt AD_mean.txt DP_mean.txt VAF.txt SF.txt GT_discordances.txt > FORMAT_JOIN.txt
		paste -d "" FORMAT_JOIN.txt FORMAT_SF.txt > FORMAT_SAMPLE.txt
		# FORMAT="GT:AD:DP:VAF:SF:GD:DV_GT:DV_DP:DV_VD:DR_GT:DR_DP:DR_VD:GK_GT:GK_DP:GK_VD"
		bcftools view -H ${sample}.merged.vcf.gz | cut -f 1-8 | sed 's/\$/\\tGT:AD:DP:VAF:SF:GD:DV_GT:DV_DP:DV_VD:DR_GT:DR_DP:DR_VD:GK_GT:GK_DP:GK_VD/' > VCF_CONTENT.txt


		bcftools view -h ${sample}.merged.vcf.gz | grep "##" > ${sample}.${assembly}.final.vcf

		echo "##FORMAT=<ID=SF,Number=1,Type=String,Description=\\"Software\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=GD,Number=1,Type=String,Description=\\"Genotype discordances. 0 same genotype and 1 different genotype\\">" >> ${sample}.${assembly}.final.vcf

		echo "##FORMAT=<ID=DV_GT,Number=1,Type=String,Description=\\"DeepVariant genotype\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=DV_DP,Number=1,Type=String,Description=\\"DeepVariant depth\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=DV_VD,Number=1,Type=String,Description=\\"DeepVariant Variant frequency\\">" >> ${sample}.${assembly}.final.vcf

		echo "##FORMAT=<ID=DR_GT,Number=1,Type=String,Description=\\"Dragen genotype\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=DR_DP,Number=1,Type=String,Description=\\"Dragen depth\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=DR_VD,Number=1,Type=String,Description=\\"Dragen Variant frequency\\">" >> ${sample}.${assembly}.final.vcf

		echo "##FORMAT=<ID=GK_GT,Number=1,Type=String,Description=\\"GATK genotype\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=GK_DP,Number=1,Type=String,Description=\\"GATK depth\\">" >> ${sample}.${assembly}.final.vcf
		echo "##FORMAT=<ID=GK_VD,Number=1,Type=String,Description=\\"GATK Variant frequency\\">" >> ${sample}.${assembly}.final.vcf

		bcftools view -h ${sample}.merged.vcf.gz | tail -n 1 | cut -f 1-10 | sed "s/.DV//" >> ${sample}.${assembly}.final.vcf

		paste -d "\\t" VCF_CONTENT.txt FORMAT_SAMPLE.txt >> ${sample}.${assembly}.final.vcf


		bgzip -c ${sample}.${assembly}.final.vcf > ${sample}.${assembly}.final.vcf.gz
		tabix -p vcf ${sample}.${assembly}.final.vcf.gz
		"""
}






process LOCALVCF {	
	label "bioinfotools"

	input:
		path inputdir
		path ref
		val sample2analyce
		
	output:
		tuple \
			val(sample2analyce_config), \
			path("${sample2analyce_config}.biallelic.vcf.gz"), emit: vcf

	script:
		sample2analyce_config = sample2analyce[0]
		"""
		if [ -f ${inputdir}/${sample2analyce_config}*.vcf.gz ]; then
			cp ${inputdir}/${sample2analyce_config}*.vcf.gz  ${sample2analyce_config}.vcf.gz
		else
			bgzip -c ${inputdir}/${sample2analyce_config}*.vcf > ${sample2analyce_config}.vcf.gz
		fi
		tabix -p vcf ${sample2analyce_config}.vcf.gz

		bcftools norm -m - -c s -f ${ref} ${sample2analyce_config}.vcf.gz | bcftools view --min-ac=1 -O z -o ${sample2analyce_config}.biallelic.vcf.gz 
		"""
}







process FORMAT2INFO {
	label "bioinfotools"

	input:
		tuple val(sample), path(final_vcf)

	output:
		tuple \
			val(sample), \
			path("${sample}.vcf_to_annotate.vcf.gz"), \
			path("${sample}.vcf_to_annotate.vcf.gz.tbi"), \
			path("${sample}.fields.txt"), emit: sample_info

	
	script:
	// GT:AD:DP:VAF:SF:GD:DV_GT:DV_DP:DV_VD:DR_GT:DR_DP:DR_VD:GK_GT:GK_DP:GK_VD
		"""
		FORMAT=(GT AD DP VAF SF GD DV_GT DV_DP DV_VD DR_GT DR_DP DR_VD GK_GT GK_DP GK_VD)

		bcftools view -h ${final_vcf} | grep "##" > ${sample}.vcf_to_annotate.vcf
		echo "##INFO=<ID=variant_id,Number=.,Type=String,Description=\\"variant identification\\">" >> ${sample}.vcf_to_annotate.vcf
		echo "##INFO=<ID=Original_pos,Number=.,Type=String,Description=\\"original position\\">" >> ${sample}.vcf_to_annotate.vcf
		for sample in \$(bcftools query -l ${final_vcf})
		do

			for field in \${FORMAT[@]}
			do
				echo "##INFO=<ID=\${sample}_\${field},Number=.,Type=String,Description=\\"\${sample} \${field}\\">" >> ${sample}.vcf_to_annotate.vcf
			done

			# echo "##INFO=<ID=\${sample}_GT,Number=.,Type=String,Description=\\"\${sample} Genotype\\">" >> ${sample}.vcf_to_annotate.vcf
			# echo "##INFO=<ID=\${sample}_AD,Number=.,Type=String,Description=\\"\${sample} Allelic depths for the ref and alt alleles in the order listed\\">" >> ${sample}.vcf_to_annotate.vcf
			# echo "##INFO=<ID=\${sample}_DP,Number=.,Type=String,Description=\\"\${sample} Approximate read depth (reads with MQ=255 or with bad mates are filtered)\\">" >> ${sample}.vcf_to_annotate.vcf
			# echo "##INFO=<ID=\${sample}_GQ,Number=.,Type=String,Description=\\"\${sample} Genotype Quality\\">" >> ${sample}.vcf_to_annotate.vcf
		done
		bcftools view -h ${final_vcf} | grep "#CHROM" | cut -f1-8 >> ${sample}.vcf_to_annotate.vcf

		newfields=""
		for field in \${FORMAT[@]}; do newfields="\$(echo "\${newfields}[;%SAMPLE\\_\${field}=%\${field}]")"; done
		bcftools query -f "variant_id=%CHROM\\_%POS\\_%REF\\_%ALT;Original_pos=%POS;\${newfields}\\n" -u ${final_vcf} | sed 's/;//2' | sed 's/,/_/g' > new_info.txt
		# bcftools query -f 'variant_id=%CHROM\\_%POS\\_%REF\\_%ALT;Original_pos=%POS;[;%SAMPLE\\_GT=%GT][;%SAMPLE\\_AD=%AD][;%SAMPLE\\_DP=%DP][;%SAMPLE\\_GQ=%GQ]\\n' -u ${final_vcf} | sed 's/;//2' | sed 's/,/_/g' > new_info.txt
		bcftools view -H ${final_vcf} | cut -f1-8 > old_info.txt
		paste -d ';' old_info.txt new_info.txt >> ${sample}.vcf_to_annotate.vcf
		
		bgzip ${sample}.vcf_to_annotate.vcf
		tabix -p vcf ${sample}.vcf_to_annotate.vcf.gz


		if bcftools view -h ${final_vcf} | grep -Fq hiConfDeNovo; then
			fields=",hiConfDeNovo,loConfDeNovo,variant_id,Original_pos"
		else
			fields=",variant_id,Original_pos"
		fi
		 
		for sample in \$(bcftools query -l ${final_vcf}); do for field in \${FORMAT[@]}; do fields="\$(echo "\${fields},\${sample}_\${field}")"; done; done
		echo \${fields} > ${sample}.fields.txt
		"""
}







process AUTOMAP {
	label "bioinfotools"
	publishDir "${params.output}/automap/", mode: 'copy'
	errorStrategy 'retry'
	
	input:
		tuple val(sample), path(final_vcf)
		val automap_assembly
		path projectDir


	output:
	tuple \
		val(sample), \
		path("*HomRegions*"), emit: roh_automap 
	
	script:
		if ( task.attempt == 1 )
			"""

			cp -R ${projectDir}/tasks/AutoMap .

			if [[ \$(bcftools query -l ${final_vcf} | wc -l) -gt 1 ]]; then
				for sample in \$(bcftools query -l ${final_vcf}); do

					bcftools view -s \${sample} -O v -o \${sample}.indv.vcf ${final_vcf}
					
					bash ./AutoMap/AutoMap_v1.2.sh \\
					--vcf \${sample}.indv.vcf \\
					--out . \\
					--genome ${automap_assembly}

					mv \${sample}/*HomRegions* .

				done
			
			else
					
				bcftools view -O v -o ${final_vcf}.vcf ${final_vcf}

				bash ./AutoMap/AutoMap_v1.2.sh \\
				--vcf ${final_vcf}.vcf \\
				--out . \\
				--genome ${automap_assembly}


				mv \$(bcftools query -l ${final_vcf})/*HomRegions* .

			fi
			"""
		else
			"""
			echo "Less than 10k variants (with quality)" > ${sample}_no_automap_HomRegions.txt
			"""
}







process VEP {
	label "vep"
	label "highcpu"
	label "highmem"
	// publishDir "${params.output}/snvs/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		path(dbscSNV)
		path(dbscSNV_tbi)
		path loFtool
		path exACpLI
		path(dbNSFP)
		path(dbNSFP_tbi)
		path maxEntScan
		path(cADD_INDELS)
		path(cADD_INDELS_tbi)
		path(cADD_SNVS)
		path(cADD_SNVS_tbi)
		path(kaviar)
		path(kaviar_tbi)
		path(cCRS_DB)
		path(cCRS_DB_tbi)
		path(dENOVO_DB)
		path(dENOVO_DB_tbi)
		path(cLINVAR)
		path(cLINVAR_tbi)
		path(gNOMADg)
		path(gNOMADg_tbi)
		path(gNOMADe)
		path(gNOMADe_tbi)
		path(gNOMADg_cov)
		path(gNOMADg_cov_tbi)
		path(gNOMADe_cov)
		path(gNOMADe_cov_tbi)
		path(cSVS)
		path(cSVS_tbi)
		path(mutScore)
		path(mutScore_tbi)
		path(mAF_FJD_COHORT)
		path(mAF_FJD_COHORT_tbi)
		path(spliceAI_SNV)
		path(spliceAI_SNV_tbi)
		path(spliceAI_INDEL)
		path(spliceAI_INDEL_tbi)
		path(REVEL)
		path(REVEL_tbi)
		path vep_cache
		path vep_plugins
		path vep_fasta
		path vep_fai
		path vep_gzi
		val vep_assembly
		tuple val(sample), path(final_vcf), path(sample_info), path(sample_info_index), path(sample_info_fields)
		val assembly
	
	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.vep.tsv"), emit: vep_tsv
	
	script:
		def dbscSNV_config = dbscSNV ? "--plugin dbscSNV,${dbscSNV} " : ''
		def loFtool_config = loFtool ? "--plugin LoFtool,${loFtool} " : ''
		def exACpLI_config = exACpLI ? "--plugin pLI,${exACpLI} " : ''
		def dbNSFP_config  = dbNSFP  ? "--plugin dbNSFP,${dbNSFP},\
LRT_pred,M-CAP_pred,MetaLR_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,\
FATHMM_pred,MetaRNN_pred,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,\
LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,\
phyloP470way_mammalian,phastCons470way_mammalian,GERP++_RS,Interpro_domain,GTEx_V8_eQTL_gene,GTEx_V8_eQTL_tissue " : ''
		def maxEntScan_config     = maxEntScan    ? "--plugin MaxEntScan,${maxEntScan} " : ''
		def cADD_config           = cADD_INDELS && cADD_SNVS ? "--plugin CADD,${cADD_INDELS},${cADD_SNVS} " : ''
		def kaviar_config         = kaviar         ? "--custom ${kaviar},kaviar,vcf,exact,0,AF,AC,AN " : ''
		def cCRS_DB_config        = cCRS_DB        ? "--custom ${cCRS_DB},gnomAD_exomes_CCR,bed,overlap,0 " : ''
		def dENOVO_DB_config      = dENOVO_DB      ? "--custom ${dENOVO_DB},denovoVariants,vcf,exact,0,SAMPLE_CT " : ''
		def cLINVAR_config        = cLINVAR        ? "--custom ${cLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNSIGCONF " : ''
		def gNOMADg_config        = gNOMADg        ? "--custom ${gNOMADg},gnomADg,vcf,exact,0,AF,AC,AN,nhomalt,grpmax,AF_grpmax,AC_grpmax,AF_nfe,AC_nfe,filt " : ''
		def gNOMADe_config        = gNOMADe        ? "--custom ${gNOMADe},gnomADe,vcf,exact,0,AF,AC,AN,nhomalt,grpmax,AF_grpmax,AC_grpmax,AF_nfe,AC_nfe,filt " : ''
		def gNOMADg_cov_config    = gNOMADg_cov    ? "--custom ${gNOMADg_cov},gnomADg_cov,vcf,overlap,0,median,perc_20x " : ''
		def gNOMADe_cov_config    = gNOMADe_cov    ? "--custom ${gNOMADe_cov},gnomADe_cov,vcf,overlap,0,median,perc_20x " : ''
		def cSVS_config           = cSVS           ? "--custom ${cSVS},CSVS,vcf,exact,0,AF,AC " : ''
		def mutScore_config       = mutScore       ? "--custom ${mutScore},Mut,vcf,exact,0,Score " : ''
		def mAF_FJD_COHORT_config = mAF_FJD_COHORT ? "--custom ${mAF_FJD_COHORT},FJD_MAF,vcf,exact,0,AF,AC,AF_DS_irdt,AC_DS_irdt,AF_P_eyeg,AC_P_eyeg " : ''
		def spliceAI_SNV_config   = spliceAI_SNV   ? "--custom ${spliceAI_SNV},SpliceAI_SNV,vcf,exact,0,SpliceAI " : ''
		def spliceAI_INDEL_config = spliceAI_INDEL ? "--custom ${spliceAI_INDEL},SpliceAI_INDEL,vcf,exact,0,SpliceAI " : ''
		def REVEL_config       	  = REVEL          ? "--custom ${REVEL},REVEL,vcf,exact,0,Score " : ''
		def sample_info_config    = sample_info    ? "--custom ${sample_info},SAMPLE,vcf,exact,0\$(cat ${sample_info_fields}) " : ''

		"""
		vep \\
		--cache --offline --dir_cache ${vep_cache} --dir_plugins ${vep_plugins} \\
		--refseq --species homo_sapiens --assembly ${vep_assembly} --force_overwrite --use_transcript_ref \\
		--verbose --fork ${params.highcpu} --tab --format vcf --no_stats \\
		--fasta ${vep_fasta} \\
		--input_file ${final_vcf} \\
		--output_file ${sample}.${assembly}.vep.tsv \\
		--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \\
		--sift p --polyphen p --allele_number --variant_class --pubmed \\
		${dbscSNV_config}\\
		${loFtool_config}\\
		${exACpLI_config}\\
		${dbNSFP_config}\\
		${maxEntScan_config}\\
		${cADD_config}\\
		${kaviar_config}\\
		${cCRS_DB_config}\\
		${dENOVO_DB_config}\\
		${cLINVAR_config}\\
		${gNOMADg_config}\\
		${gNOMADe_config}\\
		${gNOMADg_cov_config}\\
		${gNOMADe_cov_config}\\
		${cSVS_config}\\
		${mutScore_config}\\
		${mAF_FJD_COHORT_config}\\
		${spliceAI_SNV_config}\\
		${spliceAI_INDEL_config}\\
		${REVEL_config}\\
		${sample_info_config}

		"""
}


// YBQ: si el vcf es más grande de 100M, lo dividimos en trozos (de aprox 100M)

process SPLIT_VEP_TSV {
	
	input:
		tuple val(sample), path(final_vcf), path(vep_tsv)

	output:
		tuple \
			val(sample), path("split_tsv/*"), emit: split_tsv
		
	script:

		"""
		#max_size=\$((100 * 1024 * 1024)) # 100 MB in bytes
		max_size=\$((80 * 1024 * 1024)) # 80 MB in bytes
		#max_size=\$((341047/2)) 
		# Check the size of the VCF file
		vcf_size=${final_vcf.size()}

		# Calculate the number of chunks
		num_chunks=\$(( (\${vcf_size} + \${max_size} - 1) / \${max_size} ))

		if (( \${num_chunks} <= 1 )); then
			echo "The VCF file is less than 100MB. No splitting needed."
			# Prepare output directory
			output_dir="./split_tsv"
			mkdir -p \${output_dir}

			mv ${vep_tsv} "\${output_dir}/${sample}_chunk_00.tsv" 

			exit 0
		fi

		echo "VCF file is \$((\${vcf_size} / 1024 / 1024))MB. Splitting into \${num_chunks} parts."

		# Prepare output directory
		output_dir="./split_tsv"
		mkdir -p \${output_dir}

		# Extract the header from the TSV file
		header=\$(grep '^#' ${vep_tsv})

		# Calculate the number of non-header lines per chunk
		total_lines=\$(grep -vc '^#' ${vep_tsv})
		lines_per_chunk=\$(( (\${total_lines} + \${num_chunks} - 1) / \${num_chunks} + 10 ))

		echo \${lines_per_chunk}

		# Split the TSV file
		chunk=1

		grep -v '^#' ${vep_tsv} | split -l \${lines_per_chunk} - "\${output_dir}/${sample}_temp_chunk_"

		for temp_file in \${output_dir}/${sample}_temp_chunk_*; do
			current_file="\${output_dir}/${sample}_chunk-\$(printf '%03d' \${chunk}).tsv"
			echo "\${header}" > "\${current_file}"
			cat \${temp_file} >> \${current_file}
			chunk=\$((\${chunk} + 1))
			rm \${temp_file}
		done

		"""
}



process PVM {
	label "bioinfotools"
	label "highmem"
	publishDir "${params.output}/snvs/", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(vep_tsv), path(roh_automap)
		path dbNSFP_gene 
		path omim 
		path regiondict 
		path domino 
		path gene_panels
		path tissue_expression
		val maf 
		path genefilter 
		path glowgenes 
		val assembly
		path projectDir

	output:
		//tuple \
		//	val(sample), \
		//	path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"), \
		//	path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv.xlsx"), emit: pvm_tsv 
	
		tuple \
			val(sample), \
			path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"), emit: pvm_tsv 
		
		tuple \
			val(sample), \
			path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv.xlsx"), emit: pvm_xlsx

	script:
	
		def omim_field       = omim       ? "--omim ${omim} " : ''
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''
		def gene_panels_field  = gene_panels  ? "--panels ${gene_panels} " : ''

		"""
		header_row="\$(head -n 1000 ${vep_tsv} | grep "#Uploaded_variation" -n | sed 's/:.*//')"

		Rscript ${projectDir}/tasks/post-VEP_modification.R \\
		--input ${vep_tsv} \\
		--output ${sample}.${assembly}.SNV.INDEL.annotated.final.tsv \\
		--numheader \${header_row} \\
		--dbNSFPgene ${dbNSFP_gene} \\
		--regiondict ${regiondict} \\
		--domino ${domino} \\
		--expression ${tissue_expression} \\
		--automap ./ \\
		--maf ${maf} \\
		${omim_field}\\
		${genefilter_field}\\
		${glowgenes_field}\\
		${gene_panels_field}
		"""
}


process MERGE_PVM_TSV {

	publishDir "${params.output}/snvs/", mode: 'copy'
	
	input:
		tuple val(sample), path(pvm_tsvs)
		val assembly

	output:
		
		tuple \
			val(sample), path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"), emit: annotated_tsv

		
	script:
		"""
		header=\$(head -1 ${pvm_tsvs[0]})

		echo "\$header" > "${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"

		for tsv in ${pvm_tsvs}; do
			tail -n +2 \${tsv} >> "${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"
		done

		"""
}





















process BEDPROCCESING {
	label "bioinfotools"
	publishDir "${params.output}/cnvs/", mode: 'copy'
	
	input:
		path bed
		val min_target 
		val window
		val chrx
		path fai
		path projectDir


	output:
		path("*cnv.bed"), emit: bed 
	
	script:
		if ( chrx && window )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed

			python ${projectDir}/tasks/CNV_windowSize.py \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.chrx.cnv.bed

			rm \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			"""

		else if ( chrx )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed
			"""

		else if ( window )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="MT" && \$1!="chrM" && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.cnv.bed

			python ${projectDir}/tasks/CNV_windowSize.py \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.cnv.bed

			rm \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			"""

		else
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="MT" && \$1!="chrM" && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | bedtools sort -g ${fai} -i stdin > \${panel}.min${min_target}bp.cnv.bed
			"""
}







process BEDPROCCESING2 {
	label "bioinfotools"
	publishDir "${params.output}/cnvs/", mode: 'copy'
	
	input:
		path bed
		val min_target 
		val window
		val chrx
		path fai
		// path projectDir
		path ref
		path index
		path dict
		path gzi


	output:
		path("*cnv.bed"), emit: bed 
	
	script:
		if ( chrx && window )
			"""
			panel="\$(basename ${bed} .bed)"

			gatk PreprocessIntervals \
				-R ${ref} \
				-L ${bed} \
				--bin-length ${window} \
				--padding 0 \
				-imr OVERLAPPING_ONLY \
				-O \${panel}.interval_list.bed

			grep -v "@" | awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed

			// python ${projectDir}/tasks/CNV_windowSize.py \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.chrx.cnv.bed

			rm \${panel}.min${min_target}bp.chrx.cnv.bed \${panel}.min${min_target}bp.chrx.cnv.bed_unsorted
			"""

		else if ( chrx )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && (\$1=="chrX" || \$1=="X") ){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.chrx.cnv.bed
			"""

		else if ( window )
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="MT" && \$1!="chrM" && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | sort -V -k1,1 -k2,2 > \${panel}.min${min_target}bp.cnv.bed

			// python ${projectDir}/tasks/CNV_windowSize.py \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			sort -V -k1,1 -k2,2 \${panel}.min${min_target}bp.cnv.bed_unsorted | uniq > \${panel}.window125bp.min${min_target}bp.cnv.bed

			rm \${panel}.min${min_target}bp.cnv.bed \${panel}.min${min_target}bp.cnv.bed_unsorted
			"""

		else
			"""
			panel="\$(basename ${bed} .bed)"
			awk '{if((\$3-\$2)>${min_target} && \$1!="MT" && \$1!="chrM" && \$1!="chrX" && \$1!="chrY" && \$1!="Y" && \$1!="X"){print \$0}}' ${bed} | bedtools sort -g ${fai} -i stdin > \${panel}.min${min_target}bp.cnv.bed
			"""
}
















process EXOMEDEPTH {
	label "bioinfotools"
	publishDir "${params.output}/cnvs/exomedepth", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path projectDir

	output:
		tuple \
			val(runname), \
			path("exomedepth*"), emit: cnvs 

		/*tuple \
			val(runname), \
			path("exomedepth.toAnnotate.txt"), emit: toannotate, optional: true*/
	
	script:
		if ( task.attempt == 1 )
			"""
			Rscript ${projectDir}/tasks/exomeDepth.R -d . -o . -b ${bed} -n ${runname}
			"""
		else
			"""
			echo "No ExomeDepth" > exomedepth.txt
			"""
}		








process CONVADING {
	label "bioinfotools"
	publishDir "${params.output}/cnvs/convading", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path fai
		path projectDir


	output:
		tuple \
			val(runname), \
			path("CoNVaDING*"), emit: cnvs 

/*		tuple \
			val(runname), \
			path("CoNVaDING.toAnnotate.txt"), emit: toannotate, optional: true
*/	
	script:
		if ( task.attempt == 1 )
			"""
			python ${projectDir}/tasks/CoNVading_pipeline.py . ${bed} ./ ${runname} ${projectDir}/tasks ${fai}
			"""
		else
			"""
			echo "No CoNVaDING" > CoNVaDING.txt
			"""
}





process PANELCNMOPS {
	label "bioinfotools"
	publishDir "${params.output}/cnvs/panelmops", mode: 'copy'
	errorStrategy 'retry'

	input:
		path("")
		path("") 
		path bed
		val runname 
		path projectDir

	output:
		tuple \
			val(runname), \
			path("panelcn.MOPS*"), emit: cnvs 
		
		/*tuple \
			val(runname), \
			path("panelcn.MOPS.toAnnotate.txt"), emit: toannotate, optional: true*/
	
	script:
		if ( task.attempt == 1 )
			"""
			Rscript ${projectDir}/tasks/panelcnMops.R -d . -o . -b ${bed} -n ${runname}
			"""
		else
			"""
			echo "No panelcn.MOPS" > panelcn.MOPS.txt
			"""
}





process CNV_RESULT_MIXER {
	label "bioinfotools"
	publishDir "${params.output}/cnvs", mode: 'copy'
	errorStrategy 'ignore'

	input:
		/tuple val(runname), path(""), path(""), path(""), path("") /
		tuple val(runname), path(""), path(""), path("")
		path samples2analyce
		path projectDir


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.merged.bed"), emit: merged_bed 

		tuple \
			val(runname), \
			path("colnames.txt"), emit: colnames
	
	script:
		def samples2analyce_field   = samples2analyce ? "--samples ${samples2analyce} " : ""
		"""
		Rscript ${projectDir}/tasks/CNV_result_mixer2.R \
		--inputdir . \
		--outputfile ${runname}.CNV.merged.bed \
		${samples2analyce_field}

		"""
}






process ANNOTSV {
	label "annotsv"
	label "highcpu"
	label "highmem"
	// publishDir "${params.output}/cnvs", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(runname), path(merged_cnv)
		path genelist
		val annotsv_assembly
		path annotsv_path


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.annotated.tsv"), emit: annotated_cnv 
	
	script:
		// def genelist_field   = genelist ? "-candidateGenesFile ${genelist} -candidateGenesFiltering " : ""

		"""
		export ANNOTSV=${annotsv_path}

		${annotsv_path}/bin/AnnotSV \
		-SVinputFile ${merged_cnv} \
		-outputFile ${runname}.CNV.annotated.tsv \
		-annotationMode both \
		-genomeBuild ${annotsv_assembly} \
		-svtBEDcol 4 \
		-samplesidBEDcol 5 \
		-SVminSize 20 
		
		mv *_AnnotSV/${runname}.CNV.annotated.tsv .
		"""
}





process PAM {
	label "bioinfotools"
	publishDir "${params.output}/cnvs", mode: 'copy'
	errorStrategy 'ignore'
	
	input:
		tuple val(runname), path(annotated_cnv)
		tuple val(runname), path(colnames)
		path genefilter 
		path glowgenes
		path projectDir

	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.annotated.final.tsv"), emit: annotated_cnv 
	
	script:
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''

		"""
		Rscript ${projectDir}/tasks/post-AnnotSV_modification.R \\
		--input ${annotated_cnv} \\
		--outputfile ${runname}.CNV.annotated.final.tsv \\
		--extracolnames ${colnames} \\
		${genefilter_field} \\
		${glowgenes_field}
		"""
}




















process GVCF_HAPLOTYPECALLER {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	input:
		tuple val(sample), path(bam), path(bai)
		path bed
		val intervals
		val padding
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(sample), \
			path("${sample}.g.vcf"), \
			path("${sample}.g.vcf.idx"), emit: gvcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCaller" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCaller" : ""
		def intervals_field = intervals ? "-L ${bed} -ip ${padding}" : ""

		"""
		${scratch_mkdir}

		gatk --java-options "-Xmx${params.mediummem}g" \
		HaplotypeCaller ${scratch_field} \
		-R ${ref} \
		-I ${bam} \
		-ERC GVCF \
		-O ${sample}.g.vcf \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality \
		-A MappingQualityRankSumTest -A ReadPosRankSumTest \
		-A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet \
		--annotate-with-num-discovered-alleles true \
		${intervals_field}
		"""
}






process COMBINE_GVCF {	
	label "gatk"
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	input:
		path("")
		path("")
		val runname
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.g.vcf"), \
			path("${runname}.g.vcf.idx"), emit: gvcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CombineGVCFs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CombineGVCFs" : ""

		"""
		${scratch_mkdir}

		gatk CombineGVCFs ${scratch_field} \
		-R ${ref} \
		--variant *.g.vcf \
		-O ${runname}.g.vcf
		"""
}





process GENOTYPE_GVCF {	
	label "gatk"

	input:
		tuple val(runname), path(gvcf), path(idx)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.vcf"), \
			path("${runname}.vcf.idx"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_GenotypeGVCFs" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_GenotypeGVCFs" : ""
		def ped_field = ped ? "-ped ${ped}" : ""

		"""
		${scratch_mkdir}

		gatk GenotypeGVCFs ${scratch_field} \
		-R ${ref} \
		-V ${gvcf} \
		-O ${runname}.vcf \
		${ped_field}
		"""
}




process CALCULATE_GENOTYPE_POSTERIORS {	
	label "gatk"

	input:
		tuple val(runname), path(vcf), path(idx)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("filtered_INDEL_SNP_data_GP_${runname}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CalculateGenotypePosteriors" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CalculateGenotypePosteriors" : ""
		"""
		${scratch_mkdir}

		gatk CalculateGenotypePosteriors ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		-ped ${ped} \
		-O filtered_INDEL_SNP_data_GP_${runname}.vcf \
		--skip-population-priors 
		"""
}







process VARIANT_FILTRATION {	
	label "gatk"

	input:
		tuple val(runname), path(vcf)
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("filtered_INDEL_SNP_data_GP_GQfiltered_${runname}.vcf"), emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_VariantFiltration" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_VariantFiltration" : ""
		"""
		${scratch_mkdir}

		gatk VariantFiltration ${scratch_field} \
		-R ${ref} \
		-V ${vcf} \
		--filter-expression "GQ < 20" \
		--filter-name "lowGQ" \
		-O filtered_INDEL_SNP_data_GP_GQfiltered_${runname}.vcf
		"""
}







process VARIANT_ANNOTATOR {	
	label "gatk"

	input:
		tuple val(runname), path(vcf)
		path ped
		path ref
		path index
		path dict
		path reference_gzi
		path scratch
		
	output:
		tuple \
			val(runname), \
			path("${runname}.va.vcf"), 
			path("${runname}.va.vcf.idx"),emit: vcf

	script:
		def scratch_field   = scratch ? "--tmp-dir ${scratch}/${runname}_CalculateGenotypePosteriors" : ""	
		def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${runname}_CalculateGenotypePosteriors" : ""
		"""
		${scratch_mkdir}

		gatk VariantAnnotator ${scratch_field} \
		-R ${ref} \
		-V ${vcf}  \
		-O ${runname}.va.vcf \
		-A PossibleDeNovo  \
		-A StrandBiasBySample  \
		-A AS_FisherStrand \
		-ped ${ped} 
		"""
}





process MANTA {	
	label "manta"
	label "highcpu"
	label "highmem"
	publishDir "${params.output}/cnvs", mode: 'copy'

	input:
		path("")
		path("")
		path ref
		path index
		path gzi
		val runname
		
	output:
		tuple \
			val(runname), \
			path("results/variants/diploidSV.vcf.gz"), 
			path("results/variants/diploidSV.vcf.gz.tbi"), emit: diploidsv

		tuple \
			val(runname), \
			path("results/variants/candidateSV.vcf.gz"), 
			path("results/variants/candidateSV.vcf.gz.tbi"), emit: candidatesv

		tuple \
			val(runname), \
			path("results/variants/candidateSmallIndels.vcf.gz"), 
			path("results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: candidatesmallindels

	script:
		"""

		all_bams=""
		var="--bam"
		for i in *bam; do all_bams=`echo \${var} \${i} \${all_bams}`; done


		configManta.py \${all_bams} --referenceFasta ${ref} --runDir ./
		
		./runWorkflow.py -j 12

		"""
}




process ANNOTSV_VCF {
	label "annotsv"
	label "highcpu"
	label "highmem"
	publishDir "${params.output}/cnvs", mode: 'copy'
	
	input:
		tuple val(runname), path(diploidsv_vcf), path(diploidsv_index)
		val annotsv_assembly
		path annotsv_path


	output:
		tuple \
			val(runname), \
			path("${runname}.CNV.annotated.tsv"), emit: annotated_cnv 
	
	script:

		"""
		export ANNOTSV=${annotsv_path}

		${annotsv_path}/bin/AnnotSV \
		-SVinputFile ${diploidsv_vcf} \
		-outputFile ${runname}.CNV.annotated.tsv \
		-annotationMode both \
		-genomeBuild ${annotsv_assembly} \
		-SVminSize 20 
		
		mv *_AnnotSV/${runname}.CNV.annotated.tsv .
		"""
}























process PREPROCESSINTERVALS {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		path ref
		path index
		path dict
		path gzi
		path bed
		val runname
		
	output:
		tuple \
			val(runname), \
			path("${runname}.preprocessed.interval_list"), emit: interval_list

	script:
		"""
		gatk PreprocessIntervals \
			-R ${ref} \
			-L ${bed} \
			--bin-length 0 \
			-imr OVERLAPPING_ONLY \
			-O ${runname}.preprocessed.interval_list
		"""
}




process COLLECTREADCOUNTS {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(bam), path(bai)
		path ref
		path index
		path dict
		path gzi
		tuple val(runname), path(interval_list)
		
	output:
		tuple \
			val(sample), \
			path("${sample}.counts.tsv"), emit: counts

	script:
		"""
		gatk CollectReadCounts \
			-L ${interval_list} \
			-R ${ref} \
			-imr OVERLAPPING_ONLY \
			-I ${bam} \
			--format TSV \
			-O ${sample}.counts.tsv
		"""
}






process ANNOTATEINTERVALS {	
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		path ref
		path index
		path dict
		path gzi
		tuple val(runname), path(interval_list)

		
	output:
		tuple \
			val(runname), \
			path("${runname}.annotated.mytsv"), emit: annotated_intervals

	script:
		"""
		gatk AnnotateIntervals \
			-L ${interval_list} \
			-R ${ref} \
			-imr OVERLAPPING_ONLY \
			-O ${runname}.annotated.mytsv

		"""
}


process FILTERINTERVALS {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		tuple val(runname), path(interval_list)
		tuple val(runname), path(annotated_intervals)
		path("")
		
	output:
		tuple \
			val(runname), \
			path("cohort.gc.filtered.interval_list"), emit: filter_intervals

	script:
		"""
		tagI=""
		for count in *.counts.tsv
		do
			tagI="\${tagI} -I \${count}"
		done


		gatk FilterIntervals \
			-L ${interval_list} \
			--annotated-intervals ${annotated_intervals} \
			\${tagI}  \
			-imr OVERLAPPING_ONLY \
			-O cohort.gc.filtered.interval_list

		"""
}


process INTERVALLISTTOOLS {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		tuple val(runname), path(filter_intervals)
		
	output:
		path("scatter/*"), emit: scatterout

	script:
		"""
		mkdir scatter

		gatk IntervalListTools \
			--INPUT ${filter_intervals} \
			--SUBDIVISION_MODE INTERVAL_COUNT \
			--SCATTER_CONTENT 5000 \
			--OUTPUT ./scatter/
		"""
}


process DETERMINEGERMLINECONTIGPLOIDY {
	label "gatk"
	label "highcpu"
	label "highmem"
	errorStrategy 'ignore'

	input:
		tuple val(runname), path(filter_intervals)
		path("")
		path contig_ploidy_priors
		
	output:
		path("ploidy-calls"), emit: ploidy_calls
		path("ploidy-model"), emit: ploidy_model
		path("sample_index_list.txt"), emit: sample_index_list

	script:
		"""
		tagI=""
		for count in *.counts.tsv
		do
			tagI="\${tagI} -I \${count}"
		done

		gatk DetermineGermlineContigPloidy \
			-L ${filter_intervals} \
			--interval-merging-rule OVERLAPPING_ONLY \
			\${tagI} \
			--contig-ploidy-priors ${contig_ploidy_priors} \
			--output . \
			--output-prefix ploidy \
			--verbosity DEBUG

		ls ploidy-calls | sed 's/SAMPLE_//' > sample_index_list.txt
		"""
}


process GERMLINECNVCALLER {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'

	input:
		path ploidy_calls
		tuple val(runname), path(annotated_intervals)
		path("")
		path scatterdir
		
	output:
		tuple \
			val(runname), \
			path("temp*calls"),
			path("temp*model"), emit: gcnv_dir

	script:
		"""
		tagI=""
		for count in *.counts.tsv
		do
			tagI="\${tagI} -I \${count}"
		done

		prefix=\$(basename ${scatterdir})
		mkdir cohort 

		gatk GermlineCNVCaller \
			--run-mode COHORT \
			-L ${scatterdir}/scattered.interval_list \
			\${tagI} \
			--contig-ploidy-calls ${ploidy_calls} \
			--annotated-intervals ${annotated_intervals}  \
			--interval-merging-rule OVERLAPPING_ONLY \
			--output . \
			--output-prefix \${prefix} \
			--verbosity DEBUG
		"""
}


process POSTPROCESSGERMLINECNVCALLS {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'
	publishDir "${params.output}/cnvs/gatk", mode: 'copy'

	input:
		path ploidy_calls
		val sample_index
		path("")
		path("")
		path dict
		
	output:
		tuple \
			val(sample_index), \
			path("genotyped-intervals-cohort90*"), 
			path("genotyped-segments-cohort90*"),
			path("*_denoised_copy_ratios.tsv"), emit: cnv_out

	script:
		n = sample_index[0]
		"""
		tagC=""
		for  file in temp_00*-calls
		do
			tagC="\${tagC} --calls-shard-path \${file}"
		done


		tagM=""
		for  file in temp_00*-model
		do
			tagM="\${tagM} --model-shard-path \${file}"
		done


		gatk PostprocessGermlineCNVCalls \
			\${tagM} \
			\${tagC} \
			--allosomal-contig chrX --allosomal-contig chrY \
			--contig-ploidy-calls ${ploidy_calls} \
			--sample-index ${n} \
			--output-genotyped-intervals genotyped-intervals-cohort90-twelve-${n}.vcf.gz \
			--output-genotyped-segments genotyped-segments-cohort90-${n}.vcf.gz \
			--sequence-dictionary ${dict} \
			--output-denoised-copy-ratios sample_${n}_denoised_copy_ratios.tsv \
			--autosomal-ref-copy-number 2


		#bcftools query -l genotyped-intervals-cohort90-twelve-${n}.vcf.gz

		"""
}



process VCF2BED {
	label "bioinfotools"
	label "mediumcpu"
	label "mediummem"
	errorStrategy 'ignore'
	publishDir "${params.output}/cnvs/gatk", mode: 'copy'

	input:
		path("")
		val runname
		
	output:
		tuple \
			val(runname), \
			path("gatk.toAnnotate.txt"), emit: cnvs

	script:
		"""
		printf "CHR\\tSTART\\tEND\\tCNV_TYPE\\tSAMPLE\\n" > gatk.toAnnotate.txt
		for vcf in *vcf.gz
		do
		bcftools query -e 'ALT="."' -f '%CHROM\\t%POS0\\t%END\\t%ALT\\t[%SAMPLE]\\n' \${vcf} | sed 's/[<>]//g' >> gatk.toAnnotate.txt
		done
		"""
}

// ####################################################################
// ########## YBQ: Añadimos el análisis de STRs ########## 
// ####################################################################

// YBQ: necesitamos que el bai sea sí o sí .bam.bai

process EXPANSIONHUNTER {
    
    container "https://depot.galaxyproject.org/singularity/expansionhunter:5.0.0--hf366f20_0"
	publishDir "${params.output}/strs/", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path ref_fasta
	path ref_fai
	path(variant_catalog)

    output:
    tuple val(sample), path("*.vcf.gz")        , emit: vcf
    tuple val(sample), path("*.json.gz")       , emit: json
    tuple val(sample), path("*_realigned.bam") , emit: bam


    """

	bam_name=\$(basename ${bam} .bam)
	mv ${bai} \$bam_name.bam.bai

    ExpansionHunter \\
        --reads ${bam} \\
        --output-prefix ${sample} \\
        --reference ${ref_fasta} \\
        --variant-catalog ${variant_catalog}

    bgzip ${sample}.vcf
    bgzip ${sample}.json
    
	"""

}



// ####################################################################
// ########## YBQ: Añadimos el análisis del ADN mitocondrial ########## 
// ####################################################################

process PRINTREADS_CHRM {
	label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }


	input:
		tuple val(sample), path(bam), path(bai)
		path ref_fasta
		path ref_fai
		path dict
		//path gzi

	output:
		tuple \
			val(sample), \
			path("${sample}.chrM.bam"), \
			path("${sample}.chrM.bai"), emit: bam

	script:
		"""
		gatk --java-options "-Xmx${params.mediummem}g" PrintReads \
			-R ${ref_fasta} \
			-L chrM \
			--read-filter MateOnSameContigOrNoMappedMateReadFilter \
			--read-filter MateUnmappedAndUnmappedReadFilter \
			-I ${bam} \
			--read-index ${bai} \
			-O ${sample}.chrM.bam
		"""
}


process REVERTSAM {
	label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(bam), path(bai)

	output:
		tuple \
			val(sample), \
			path("${sample}.chrM.unmapped.bam"), emit: bam

	script:
		"""
		gatk --java-options "-Xmx${params.mediummem}g" RevertSam \
			INPUT=${bam} \
			OUTPUT=${sample}.chrM.unmapped.bam \
			OUTPUT_BY_READGROUP=false \
			VALIDATION_STRINGENCY=LENIENT \
			ATTRIBUTE_TO_CLEAR=FT \
			ATTRIBUTE_TO_CLEAR=CO \
			SORT_ORDER=queryname \
			RESTORE_ORIGINAL_QUALITIES=false
		"""
}


process SAMTOFASTQ {
	label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(bam)

	output:
		tuple \
			val(sample), \
			path("${sample}_R1.chrM.fastq"), \
			path("${sample}_R2.chrM.fastq"), emit: fastq

	script:
		"""
		gatk --java-options "-Xmx${params.mediummem}g" SamToFastq \
			INPUT=${bam} \
			FASTQ=${sample}_R1.chrM.fastq \
			SECOND_END_FASTQ=${sample}_R2.chrM.fastq \
			NON_PF=true
		"""
}



process BWA_CHRM {
	label "bwa"
	label "highcpu"
	label "highmem"

	input:
		tuple val(sample), path(forward), path(reverse)
		val ref_name
		path index
		path ref

	output:
		tuple \
			val(sample), \
			path("${sample}.${ref_name}.mapped.bam"), emit: mapped_bam

	script:
		// sample  = fastq[0]
		// forward = fastq[1]
		// reverse = fastq[2]
		"""
		bwa mem \\
		${ref} \\
		${forward} \\
		${reverse} > ${sample}.${ref_name}.mapped.bam

		"""
}


process MERGEBAMALIGNMENT_CHRM {
	label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(mapped_bam), path(unmapped_bam)
		path ref
		path index
		path dict
		path scratch
		val ref_name

	output:
		tuple \
			val(sample), \
			path("${sample}.${ref_name}.mapped.merged.bam"), emit: merged_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}_MergeBamAlignment" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}_MergeBamAlignment" : ""

		"""

		${scratch_mkdir}

		gatk --java-options "-Xmx${params.mediummem}g" MergeBamAlignment ${scratch_field} \
			--VALIDATION_STRINGENCY SILENT \
			--EXPECTED_ORIENTATIONS FR \
			--ATTRIBUTES_TO_RETAIN X0 \
			--ATTRIBUTES_TO_REMOVE NM \
			--ATTRIBUTES_TO_REMOVE MD \
			--ALIGNED_BAM ${mapped_bam} \
			--UNMAPPED_BAM ${unmapped_bam} \
			--OUTPUT ${sample}.${ref_name}.mapped.merged.bam \
			--REFERENCE_SEQUENCE ${ref} \
			--PAIRED_RUN true \
			--SORT_ORDER "unsorted" \
			--IS_BISULFITE_SEQUENCE false \
			--ALIGNED_READS_ONLY false \
			--CLIP_ADAPTERS false \
			--MAX_RECORDS_IN_RAM 2000000 \
			--ADD_MATE_CIGAR true \
			--MAX_INSERTIONS_OR_DELETIONS -1 \
			--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
			--UNMAPPED_READ_STRATEGY COPY_TO_TAG \
			--ALIGNER_PROPER_PAIR_FLAGS true \
			--UNMAP_CONTAMINANT_READS true \
			--ADD_PG_TAG_TO_READS false
		"""
}


process MARKDUPLICATES_CHRM {
	label "gatk"
	label "mediumcpu"
	label "mediummem"
	maxRetries 3

	input:
		tuple val(sample), path(merged_bam)
		path scratch
		val ref_name

	output:
		tuple \
			val(sample), \
			path("${sample}.${ref_name}.dedupped.sorted.bam"), emit: deduppedsorted_bam
		tuple \
			val(sample), \
			path("marked_dup_metrics_${sample}.${ref_name}.txt"), emit: dedupped_txt

	script:
		//def scratch_field = scratch ? "--conf 'spark.local.dir=${scratch}/${sample}.${ref_name}_MarkDuplicatesSpark'" : ""
		//def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}.${ref_name}_MarkDuplicatesSpark" : ""

		"""
		gatk MarkDuplicates \
		-I ${merged_bam} \
		-O ${sample}.${ref_name}.dedupped.sorted.bam \
		-M marked_dup_metrics_${sample}.${ref_name}.txt \
		--VALIDATION_STRINGENCY "SILENT" \
		--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
		--ASSUME_SORT_ORDER queryname \
		--CLEAR_DT "false" \
		--ADD_PG_TAG_TO_READS false 
		

		#chmod 777 \$(find . -user root) 
		chmod 777 *.dedupped.sorted.bam* .command.trace marked_dup_metrics*
		"""
}

process SORTSAM_CHRM {
	label "gatk"	
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(merged_bam)
		path scratch
		val ref_name
		
	output:
		tuple \
			val(sample), \
			path("${sample}.${ref_name}.sorted.bam"), \
			path("${sample}.${ref_name}.sorted.bai"), emit: sorted_bam

	script:
		def scratch_field = scratch ? "--TMP_DIR ${scratch}/${sample}.${ref_name}_SortSam" : ""
		def scratch_mkdir = scratch ? "mkdir -p ${scratch}/${sample}.${ref_name}_SortSam" : ""
		"""

		${scratch_mkdir}

		gatk SortSam ${scratch_field} \
		--INPUT ${merged_bam} \
		--OUTPUT ${sample}.${ref_name}.sorted.bam \
		--SORT_ORDER "coordinate" \
		--CREATE_INDEX true \
		--CREATE_MD5_FILE false 
		"""
}


process COLLECTWGSMETRICTS {
	label "gatk"	
	maxRetries 3
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
	publishDir "${params.output}/mitochondria/metrics", mode: 'copy'

	input:
		tuple val(sample), path(bam), path(bai)
		path ref
		
	output:
		tuple \
			val(sample), \
			path("metrics.txt"), \
			path("theoretical_sensitivity.txt"), emit: metrics

	script:

		"""
		gatk CollectWgsMetrics \
			--INPUT ${bam} \
			--VALIDATION_STRINGENCY "SILENT" \
			--REFERENCE_SEQUENCE ${ref} \
			--OUTPUT metrics.txt \
			--USE_FAST_ALGORITHM true \
			--COVERAGE_CAP 100000 \
			--INCLUDE_BQ_HISTOGRAM true \
			--THEORETICAL_SENSITIVITY_OUTPUT theoretical_sensitivity.txt
		"""
}

process GATK4_MUTECT2 {

    label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(bam), path(bai)
	path ref_fasta
	path ref_fai
	path dict
	path scratch
	val ref_name
	val region


    output:

    tuple val(sample), path("*.vcf")     , emit: vcf
    tuple val(sample), path("*.stats")      , emit: stats


	script:
	def scratch_field   = scratch ? "--tmp-dir ${scratch}/${sample}_HaplotypeCaller" : ""	
	def scratch_mkdir   = scratch ? "mkdir -p ${scratch}/${sample}_HaplotypeCaller" : ""

	"""
	${scratch_mkdir}

	gatk --java-options "-Xmx${params.mediummem}g" \\
        Mutect2 ${scratch_field} \\
        -I ${bam} \\
		-R ${ref_fasta} \\
		--read-filter MateOnSameContigOrNoMappedMateReadFilter \\
        --read-filter MateUnmappedAndUnmappedReadFilter \\
		-O ${sample}.${ref_name}.vcf \\
		--bam-output bamout.${ref_name}.bam  \\
		-L ${region} \\
		--annotation StrandBiasBySample \\
        --mitochondria-mode \\
        --max-mnp-distance 0 \\
		--max-reads-per-alignment-start 75

    """
}


process LIFTOVER_CHRM {
	label "gatk"	
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(vcf)
		path ref_fasta
		path dict
		path chain

	output:
		tuple \
			val(sample), \
			path("*.shifted_back.vcf"), emit: shiftedback_vcf
		tuple \
			val(sample), \
			path("*.rejected.vcf"), emit: rejected_vcf

	script:

		"""
		gatk LiftoverVcf \
			--INPUT ${vcf} \
			--OUTPUT ${sample}.shifted_back.vcf \
			-R ${ref_fasta} \
			--CHAIN ${chain} \
			--REJECT ${sample}.rejected.vcf
		"""
}

process MERGE_VCFS_CHRM {
	label "gatk"	
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
		tuple val(sample), path(vcf_chrmt), path(vcf_shiftedback)

	output:
		tuple \
			val(sample), \
			path("${sample}.merged.vcf"), 
			path("${sample}.merged.vcf.idx"), emit: merged_vcf

	script:

		"""
		gatk MergeVcfs \
			--I ${vcf_chrmt} \
			--I ${vcf_shiftedback} \
			--O ${sample}.merged.vcf
		"""
}


process MERGE_MUTECT_STATS {
    label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(shifted_stats), path(non_shifted_stats)

    output:

    tuple val(sample), path("${sample}.combined.stats")   , emit: combined_stats

	script:

		"""

		gatk MergeMutectStats \\
		--stats ${shifted_stats} \\
		--stats ${non_shifted_stats} \\
		-O ${sample}.combined.stats

		"""
	
}


process FILTER_MUTECT_CALLS_INITIAL {
    label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(vcf), path(vcf_index), path(stats)
	path ref_fasta
	path ref_fai
	path dict
	path (blacklisted_sites)
	path (blacklisted_sites_index)

    output:

    tuple val(sample), path("${sample}.filtered.vcf"), path("${sample}.filtered.vcf.idx"), emit: filtered_vcf
	
	script:


		"""

		gatk --java-options "-Xmx${params.mediummem}g" \\
			FilterMutectCalls  \\
			-V ${vcf} \\
			-R ${ref_fasta} \\
			-O ${sample}.mutect_filtered.vcf \\
			--stats ${stats} \\
			--mitochondria-mode  \\
			--max-alt-allele-count 4 \\
			--min-allele-fraction 0 \\

		gatk --java-options "-Xmx${params.mediummem}g" VariantFiltration  \\
			-V 	${sample}.mutect_filtered.vcf \\
			-O  ${sample}.filtered.vcf \\
			--apply-allele-specific-filters \\
			--mask ${blacklisted_sites} \\
			--mask-name "blacklisted_site"

		"""	
}


process SPLITMULTIALLELICS_AND_REMOVENONPASS_SITES {
    label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(vcf), path(vcf_index)
	path ref_fasta
	path ref_fai
	path dict


    output:
    tuple val(sample), path("${sample}.filtered.splitAndPassOnly.vcf"), path("${sample}.filtered.splitAndPassOnly.vcf.idx")     , emit: split_pass_vcf

	script:

		"""
		gatk --java-options "-Xmx${params.mediummem}g" LeftAlignAndTrimVariants  \\
			-V ${vcf} \\
			-R ${ref_fasta} \\
			-O ${sample}.filtered.split.vcf \\
			--split-multi-allelics \\
			--dont-trim-alleles \\
			--keep-original-ac 

		gatk --java-options "-Xmx${params.mediummem}g" SelectVariants \
			-V ${sample}.filtered.split.vcf \\
			-O ${sample}.filtered.splitAndPassOnly.vcf \\
			--exclude-filtered

		"""	
}


process HAPLOCHECK {
	container "https://depot.galaxyproject.org/singularity/haplocheck:1.3.3--h4a94de4_0"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

	input:
	tuple val(sample), path(vcf), path(vcf_index)

    output:

	tuple val(sample), \
		path("${sample}.contamination.txt"), \
		path("${sample}.mean_het_major.txt"), \
		path("${sample}.mean_het_minor.txt"), emit: haplocheck_contamination

    //tuple val(sample), path("${sample}.contamination.txt") , emit: hasContamination
    //tuple val(sample), path("${sample}.major_hg.txt") , emit: major_hg
    //tuple val(sample), path("${sample}.minor_hg.txt") , emit: minor_hg
    //tuple val(sample), path("${sample}.mean_het_major.txt") , emit: major_level
    //tuple val(sample), path("${sample}.mean_het_minor.txt") , emit: minor_level
    //tuple val(sample), path("${sample}.html"), emit: html

	script:

		"""
		haplocheck --raw --out ${sample} ${vcf}
		sed 's/\"//g' ${sample}.raw.txt > ${sample}-noquotes
		grep "Sample" ${sample}-noquotes > headers

		FORMAT_ERROR="Bad contamination file format"

		if [ `awk 'BEGIN { FS = "\t" };{print \$2}' headers` != "Contamination Status" ]; then
			echo \$FORMAT_ERROR; exit 1
		fi
		if [ `awk 'BEGIN { FS = "\t" };{print \$10}' headers` != "Major Haplogroup" ]; then
			echo \$FORMAT_ERROR; exit 1
		fi
		if [ `awk 'BEGIN { FS = "\t" };{print \$12}' headers` != "Minor Haplogroup" ]; then
			echo \$FORMAT_ERROR; exit 1
		fi
		if [ `awk 'BEGIN { FS = "\t" };{print \$8}' headers` != "Major Heteroplasmy Level" ]; then
			echo \$FORMAT_ERROR; exit 1
		fi
		if [ `awk 'BEGIN { FS = "\t" };{print \$9}' headers` != "Minor Heteroplasmy Level" ]; then
			echo \$FORMAT_ERROR; exit 1
		fi

		grep -v "Sample" ${sample}-noquotes >  ${sample}.output-data
		awk -F "\t" '{print \$2}' ${sample}.output-data > ${sample}.contamination.txt
		awk -F "\t" '{print \$10}' ${sample}.output-data > ${sample}.major_hg.txt
		awk -F "\t" '{print \$12}' ${sample}.output-data > ${sample}.minor_hg.txt
		awk -F "\t" '{print \$8}' ${sample}.output-data > ${sample}.mean_het_major.txt
		awk -F "\t" '{print \$9}' ${sample}.output-data > ${sample}.mean_het_minor.txt


		"""	

}

process FILTER_MUTECT_CALLS_CONTAMINATION {
    label "gatk"
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(vcf), path(vcf_index), path(stats), path (hasContamination), path (major_level), path (minor_level)
	path ref_fasta
	path ref_fai
	path dict
	path (blacklisted_sites)
	path (blacklisted_sites_index)


    output:

    tuple val(sample), path("${sample}.contamination_filtered.vcf"), path("${sample}.contamination_filtered.vcf.idx"), emit: contamination_filtered_vcf
	
	script:

		"""

		hasContamination_check=\$(cat ${hasContamination} )
		contamination_major=\$(cat ${major_level} )
		contamination_minor=\$(cat ${minor_level} )

		echo \$hasContamination_check \$contamination_major  \$contamination_minor 

		if [ \$hasContamination_check == "YES" ] 
		then 
			if [ \$contamination_major == "0.0" ]
			then
				hc_contamination=\$contamination_minor
			else 
				hc_contamination=\$(echo "1.0 - \$contamination_major" | bc)
			fi
		else 
			hc_contamination=0
		fi

		echo "#######"
		echo \$hc_contamination
		echo "#######"

		gatk --java-options "-Xmx${params.mediummem}g" \\
			FilterMutectCalls  \\
			-V ${vcf} \\
			-R ${ref_fasta} \\
			-O ${sample}.mutect_filtered.vcf \\
			--stats ${stats} \\
			--mitochondria-mode  \\
			--max-alt-allele-count 4 \\
			--min-allele-fraction 0 \\
			--contamination-estimate \$hc_contamination

		gatk --java-options "-Xmx${params.mediummem}g" VariantFiltration  \\
			-V 	${sample}.mutect_filtered.vcf \\
			-O  ${sample}.contamination_filtered.vcf \\
			--apply-allele-specific-filters \\
			--mask ${blacklisted_sites} \\
			--mask-name "blacklisted_site"

		"""	
}

process SPLITMULTIALLELICSSITES_CHR {
    label "gatk"
	publishDir "${params.output}/mitochondria/snvs", mode: 'copy'
	errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
	tuple val(sample), path(vcf), path(vcf_index)
	path ref_fasta
	path ref_fai
	path dict


    output:
    tuple val(sample), path("${sample}.chrM.final.vcf"), path("${sample}.chrM.final.vcf.idx")     , emit: final_chrM_vcf

	script:

		"""
		gatk --java-options "-Xmx${params.mediummem}g" LeftAlignAndTrimVariants  \\
			-V ${vcf} \\
			-R ${ref_fasta} \\
			-O ${sample}.chrM.final.vcf \\
			--split-multi-allelics \\
			--dont-trim-alleles \\
			--keep-original-ac 

		"""	
}



process FORMAT2INFO_CHRM {
	label "bioinfotools"

	input:
		tuple val(sample), path(final_vcf), path(final_vcf_index)

	output:
		tuple \
			val(sample), \
			path("${sample}.vcf_to_annotate.vcf.gz"), \
			path("${sample}.vcf_to_annotate.vcf.gz.tbi"), emit: sample_info

	
	script:

		"""
		awk 'BEGIN {FS="\t"; OFS="\t"} {
			if (\$0 ~ /^##/) {
					print \$0
				}    
				else if (\$0 ~ /^#CHROM/) {
					# Add the new ##INFO=<ID=filter> line before the #CHROM line
					print "##INFO=<ID=AD_REF,Number=1,Type=String,Description=\"AD ref\">"
					print "##INFO=<ID=AD_ALT,Number=1,Type=String,Description=\"AD alt\">"
					print "##INFO=<ID=filter,Number=.,Type=String,Description=\"filter flag\">"
					print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8
				} else {
				# Extract fields
				chrom = \$1
				pos = \$2
				id = \$3
				ref = \$4
				alt = \$5
				qual = \$6
				filter = \$7
				info = \$8
				format = \$9
				sample = \$10

				# Split FORMAT and sample fields
				split(format, format_fields, ":")
				split(sample, sample_fields, ":")

				# Create a dictionary to map FORMAT fields to sample values
				for (i in format_fields) {
					format_dict[format_fields[i]] = sample_fields[i]
				}

				# Extract desired fields
				ad = format_dict["AD"]
				af = format_dict["AF"]


				# Split AD into AD_REF and AD_ALT
				split(ad, ad_values, ",")
				ad_ref = ad_values[1]
				ad_alt = ad_values[2]

				# Append new fields to INFO
				new_info = info ";AD_REF=" ad_ref ";AD_ALT=" ad_alt ";AF=" af ";filter=" filter

				# Print the transformed line
				print chrom, pos, id, ref, alt, qual, filter, new_info
			}
		}' ${final_vcf} > ${sample}.vcf_to_annotate.vcf
		
		bgzip ${sample}.vcf_to_annotate.vcf
		tabix -p vcf ${sample}.vcf_to_annotate.vcf.gz

		"""
}






process VEP_CHRM {
	label "vep"
	label "highcpu"
	label "highmem"
	// publishDir "${params.output}/snvs/", mode: 'copy'
	//publishDir "${params.output}/mitochondria/snvs", mode: 'copy'

	input:
		path(dbscSNV)
		path(dbscSNV_tbi)
		path loFtool
		path exACpLI
		path(dbNSFP)
		path(dbNSFP_tbi)
		path maxEntScan
		path(cADD_INDELS)
		path(cADD_INDELS_tbi)
		path(cADD_SNVS)
		path(cADD_SNVS_tbi)
		path(kaviar)
		path(kaviar_tbi)
		path(cCRS_DB)
		path(cCRS_DB_tbi)
		path(dENOVO_DB)
		path(dENOVO_DB_tbi)
		path(cLINVAR)
		path(cLINVAR_tbi)
		path(cSVS)
		path(cSVS_tbi)
		path(mutScore)
		path(mutScore_tbi)
		path(mAF_FJD_COHORT)
		path(mAF_FJD_COHORT_tbi)
		path(spliceAI_SNV)
		path(spliceAI_SNV_tbi)
		path(spliceAI_INDEL)
		path(spliceAI_INDEL_tbi)
		path(REVEL)
		path(REVEL_tbi)
		path vep_cache
		path vep_plugins
		path vep_fasta
		path vep_fai
		path vep_gzi
		val vep_assembly
		tuple val(sample), path(vcf_to_annotate), path(vcf_to_annotate_idx), path(vcf_w_sample_info), path(vcf_w_sample_info_index)
		val assembly
		path chrM_polymorphisms
		path chrM_polymorphisms_tbi
		path chrM_disease
		path chrM_disease_tbi
		path mitomap_genomeloci
		path mitomap_genomeloci_tbi
		path mitotip
		path mitotip_tbi
		path gnomad_chrM_reduced
		path gnomad_chrM_reduced_tbi
		path gnomad_chrM_full
		path gnomad_chrM_full_tbi
		path mitimpact
		path mitimpact_tbi
		path t_apogee
		path t_apogee_tbi
	
	
	output:
		tuple \
			val(sample), \
			path("${sample}.${assembly}.vep.tsv"), emit: vep_tsv
	
	script:
		def dbscSNV_config = dbscSNV ? "--plugin dbscSNV,${dbscSNV} " : ''
		def loFtool_config = loFtool ? "--plugin LoFtool,${loFtool} " : ''
		def exACpLI_config = exACpLI ? "--plugin pLI,${exACpLI} " : ''
		def dbNSFP_config  = dbNSFP  ? "--plugin dbNSFP,${dbNSFP},\
LRT_pred,M-CAP_pred,MetaLR_pred,MetaSVM_pred,MutationAssessor_pred,MutationTaster_pred,PROVEAN_pred,\
FATHMM_pred,MetaRNN_pred,PrimateAI_pred,DEOGEN2_pred,BayesDel_addAF_pred,BayesDel_noAF_pred,ClinPred_pred,\
LIST-S2_pred,Aloft_pred,fathmm-MKL_coding_pred,fathmm-XF_coding_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,\
phyloP470way_mammalian,phastCons470way_mammalian,GERP++_RS,Interpro_domain,GTEx_V8_eQTL_gene,GTEx_V8_eQTL_tissue " : ''
		def maxEntScan_config     = maxEntScan    ? "--plugin MaxEntScan,${maxEntScan} " : ''
		def cADD_config           = cADD_INDELS && cADD_SNVS ? "--plugin CADD,${cADD_INDELS},${cADD_SNVS} " : ''
		def kaviar_config         = kaviar         ? "--custom ${kaviar},kaviar,vcf,exact,0,AF,AC,AN " : ''
		def cCRS_DB_config        = cCRS_DB        ? "--custom ${cCRS_DB},gnomAD_exomes_CCR,bed,overlap,0 " : ''
		def dENOVO_DB_config      = dENOVO_DB      ? "--custom ${dENOVO_DB},denovoVariants,vcf,exact,0,SAMPLE_CT " : ''
		def cLINVAR_config        = cLINVAR        ? "--custom ${cLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNSIGCONF " : ''
		def cSVS_config           = cSVS           ? "--custom ${cSVS},CSVS,vcf,exact,0,AF,AC " : ''
		def mutScore_config       = mutScore       ? "--custom ${mutScore},Mut,vcf,exact,0,Score " : ''
		def mAF_FJD_COHORT_config = mAF_FJD_COHORT ? "--custom ${mAF_FJD_COHORT},FJD_MAF,vcf,exact,0,AF,AC,AF_DS_irdt,AC_DS_irdt,AF_P_eyeg,AC_P_eyeg " : ''
		def spliceAI_SNV_config   = spliceAI_SNV   ? "--custom ${spliceAI_SNV},SpliceAI_SNV,vcf,exact,0,SpliceAI " : ''
		def spliceAI_INDEL_config = spliceAI_INDEL ? "--custom ${spliceAI_INDEL},SpliceAI_INDEL,vcf,exact,0,SpliceAI " : ''
		def REVEL_config       	  = REVEL          ? "--custom ${REVEL},REVEL,vcf,exact,0,Score " : ''
		def chrM_polymorphisms_config    = chrM_polymorphisms ?  "--custom ${chrM_polymorphisms},Mitomap_GenBank,vcf,exact,0,AC,AF " : ''
		def chrM_disease_config    = chrM_disease ?  "--custom ${chrM_disease},Mitomap_disease,vcf,exact,0,aachange,Disease,DiseaseStatus,HGFL " : ''
		def chrM_mitomapgenomeloci_config = mitomap_genomeloci ? "--custom ${mitomap_genomeloci},Mitomap_genomeloci,bed " : ''
		def chrM_mitotip_config    = mitotip ?  "--custom ${mitotip},MitoTip,vcf,exact,0,MitoTIP_Score,MitoTIP_Quartile,MitoTIP_Count,MitoTIP_Percentage,MitoTIP_Status " : ''
		def chrM_gnomad_reduced_config    = gnomad_chrM_reduced ?  "--custom ${gnomad_chrM_reduced},gnomAD,vcf,exact,0,AC_hom,AC_het,AF_hom,AF_het,AN,max_observed_heteroplasmy,FILTER " : ''
		def chrM_gnomad_full_config    = gnomad_chrM_full ?  "--custom ${gnomad_chrM_full},gnomAD_full,vcf,exact,0,hap_defining_variant " : ''
		def sample_info_config    = vcf_w_sample_info    ? "--custom ${vcf_w_sample_info},SAMPLE,vcf,exact,0,AD_REF,AD_ALT,AF,DP,filter " : ''
		def mitimpact_config    = mitimpact    ? "--custom ${mitimpact},MitImpact,vcf,exact,0,HGVS,\
Functional_effect_general,Functional_effect_detailed,PolyPhen2,SIFT,SIFT4G,VEST_pvalue,VEST,VEST_FDR,Mitoclass1,SNPDryad,FATHMM,\
AlphaMissense,CADD_score,CADD_phred_score,CADD,PROVEAN,MutationAssessor,EFIN_SP,EFIN_HD,MLC,PANTHER,CAROL,Condel,\
APOGEE1_score,APOGEE1,APOGEE2_score,APOGEE2,PhastCons_100V,PhyloP_100V,PhyloP_470Way,PhastCons_470Way,\
HelixMTdb_AC_hom,HelixMTdb_AF_hom,HelixMTdb_AC_het,HelixMTdb_AF_het,HelixMTdb_mean_ARF,HelixMTdb_max_ARF,ToMMo_54KJPN_AC,ToMMo_54KJPN_AF,ToMMo_54KJPN_AN " : ''
		def t_apogee_config    = t_apogee ?  "--custom ${t_apogee},tRNA_APOGEE,vcf,exact,0,t-APOGEE_unbiased_score " : ''


		"""
		vep \\
		--cache --offline --dir_cache ${vep_cache} --dir_plugins ${vep_plugins} \\
		--refseq --species homo_sapiens --assembly ${vep_assembly} --force_overwrite --use_transcript_ref \\
		--distance 0 --verbose --fork ${params.highcpu} --tab --format vcf --no_stats \\
		--fasta ${vep_fasta} \\
		--input_file ${vcf_to_annotate} \\
		--output_file ${sample}.${assembly}.vep.tsv \\
		--appris --canonical --ccds \\
		--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \\
		--sift p --allele_number --variant_class --pubmed \\
		--polyphen p --tsl --uniprot \\
		${dbscSNV_config}\\
		${loFtool_config}\\
		${exACpLI_config}\\
		${dbNSFP_config}\\
		${maxEntScan_config}\\
		${cADD_config}\\
		${kaviar_config}\\
		${cCRS_DB_config}\\
		${dENOVO_DB_config}\\
		${cLINVAR_config}\\
		${cSVS_config}\\
		${mutScore_config}\\
		${mAF_FJD_COHORT_config}\\
		${spliceAI_SNV_config}\\
		${spliceAI_INDEL_config}\\
		${REVEL_config}\\
		${chrM_polymorphisms_config} \\
		${chrM_mitomapgenomeloci_config} \\
		${chrM_disease_config} \\
		${chrM_mitotip_config} \\
		${chrM_gnomad_full_config} \\
		${chrM_gnomad_reduced_config} \\
		${mitimpact_config} \\
		${t_apogee_config} \\
		${sample_info_config}

		"""
}



process PVM_CHRM {
	label "bioinfotools_MT"
	label "highmem"
	publishDir "${params.output}/mitochondria/snvs", mode: 'copy'
	errorStrategy 'ignore'

	input:
		tuple val(sample), path(vep_tsv)
		path dbNSFP_gene 
		path omim 
		path regiondict 
		path domino 
		path gene_panels
		path tissue_expression
		path genefilter 
		path glowgenes 
		val assembly
		path projectDir

	output:
		//tuple \
		//	val(sample), \
		//	path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv"), \
		//	path("${sample}.${assembly}.SNV.INDEL.annotated.final.tsv.xlsx"), emit: pvm_tsv 
	
		tuple \
			val(sample), \
			path("${sample}.${assembly}.CHRM.SNV.INDEL.annotated.final.tsv"), emit: pvm_tsv 
		
		tuple \
			val(sample), \
			path("${sample}.${assembly}.CHRM.SNV.INDEL.annotated.final.tsv.xlsx"), emit: pvm_xlsx

	script:
	
		def omim_field       = omim       ? "--omim ${omim} " : ''
		def genefilter_field = genefilter ? "--genefilter ${genefilter} " : ''
		def glowgenes_field  = glowgenes  ? "--glowgenes ${glowgenes} " : ''
		def gene_panels_field  = gene_panels  ? "--panels ${gene_panels} " : ''

		"""
		header_row="\$(head -n 1000 ${vep_tsv} | grep "#Uploaded_variation" -n | sed 's/:.*//')"

		Rscript ${projectDir}/tasks/post-VEP_CHRM_modification.R \\
		--input ${vep_tsv} \\
		--output ${sample}.${assembly}.CHRM.SNV.INDEL.annotated.final.tsv \\
		--numheader \${header_row} \\
		--dbNSFPgene ${dbNSFP_gene} \\
		--regiondict ${regiondict} \\
		--domino ${domino} \\
		--expression ${tissue_expression} \\
		--automap ./ \\
		${omim_field}\\
		${genefilter_field}\\
		${glowgenes_field}\\
		${gene_panels_field}
		"""
}


// ####################################################################





