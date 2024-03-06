from subprocess import call 
from glob import glob
import sys
import argparse
import os
import numpy as np

#Arguments needed to run the script
bam_folder = sys.argv[1] + '*.bam'
bedfile = sys.argv[2]
output_folder = sys.argv[3]
fai = sys.argv[4]

#Creating output folder if doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#Creating raw_coverages folder if doesn't exist
raw_coverages_folder = output_folder + 'raw_coverages_folder'
if not os.path.exists(raw_coverages_folder):
    os.makedirs(raw_coverages_folder)

#Creating StartWithMatchScore folder if doesn't exist
StartWithMatchScore = output_folder + 'StartWithMatchScore'
if not os.path.exists(StartWithMatchScore):
    os.makedirs(StartWithMatchScore)

#Defining function in order to obtain sample name from the bam file path
def get_sample_name(sample_path):
	sample_name = sample_path[sample_path.rfind('/') +1:]
	sample_name = sample_name[:sample_name.find('_')]
	return sample_name

#Defining function to calculate the gene list of the sample
def get_gene_list(cov_file):
	gene_list = list()
	for line in cov_file:
		line = line.strip('\n').split('\t')
		if line[3].strip() not in gene_list:
			gene_list.append(line[3].strip())
	return gene_list

#Defining function to calculate the total average coverage of the sample
def get_avg_total_cov(cov_file):
	all_mean_coverages = list()
	for line in cov_file:
		line = line.strip('\n').split('\t')
		all_mean_coverages.append(float(line[4].strip()))

	sample_all_mean_cov = np.mean(all_mean_coverages)
	return sample_all_mean_cov

#Defining function to calculate the autosomal average coverage of the sample
def get_avg_auto_cov(cov_file):
	autosomal_mean_coverages = list()
	for line in cov_file:
		line = line.strip('\n').split('\t')

		if line[0].strip() != 'chrX':
			autosomal_mean_coverages.append(float(line[4].strip()))

	sample_auto_mean_cov = np.mean(autosomal_mean_coverages)
	return sample_auto_mean_cov

#Defining function to calculate the average gene coverages of the sample 
def get_avg_gene_cov(cov_file,gene_list):
	avg_gene_cov_dict = dict()
	for gene in gene_list:
		avg_gene_coverages = list()

		for line in cov_file:
			line = line.strip('\n').split('\t')

			if gene == line[3].strip():
				avg_gene_coverages.append(float(line[4].strip()))

		avg_gene_cov_dict.update({gene:np.mean(avg_gene_coverages)})
	return avg_gene_cov_dict


#Calculating mean coverage of each target region (bedfile) for every sample and storing it in the working folder
for bam_file in glob(bam_folder):
	sample_name = get_sample_name(bam_file)
	print '-----'
	print 'BedFile: ' + bedfile + '\nBamFile: ' + bam_file + '\nFai: ' + fai
	print 'Calculating coverage of sample: ' + sample_name + ' ...'
	call('bedtools coverage -a ' + bedfile + ' -b ' + bam_file + ' -sorted -g ' + fai + ' -mean| sort -k1,1 -k2,2n  > ' + raw_coverages_folder + '/' + sample_name + '_mean_coverage.txt',shell = True)
#	call('bedtools coverage -a ' + bedfile + ' -b ' + bam_file + ' -sorted -g /mnt/genetica/ionut/GeneticaPipeDB/pipeline/CoNVaDING-1.2.1/rCRS_fasta.fai -mean| sort -k1,1 -k2,2n  > ' + raw_coverages_folder + '/' + sample_name + '_mean_coverage.txt',shell = True)


	#call('bedtools coverage -a ' + bedfile + ' -b ' + bam_file + ' -mean| sort -k1,1 -k2,2n > ' + raw_coverages_folder + '/' + sample_name + '_mean_coverage.txt',shell = True)

	print 'Saved in ' + raw_coverages_folder

#Iterate over all the coverage files to append the different mean coverages
ext_cov_file_header = ['CHR','START','STOP','GENE','REGION_COV','AVG_AUTOSOMAL_COV','AVG_TOTAL_COV','AVG_GENE_COV','NORMALIZED_AUTOSOMAL','NORMALIZED_TOTAL','NORMALIZED_GENE']
for coverage_file_path in glob(raw_coverages_folder + '/*'):
	coverage_file = open(coverage_file_path, 'r')
	coverage_file = coverage_file.readlines()

	#Getting the total average coverage 
	total_mean_cov = get_avg_total_cov(coverage_file)

	#Getting the autosomal average coverage
	auto_mean_cov = get_avg_auto_cov(coverage_file)

	#Getting the average gene coverages
	gene_list = get_gene_list(coverage_file)
	gene_mean_cov = get_avg_gene_cov(coverage_file,gene_list)

	#Append the different coverages for each line and listing all lines for saving
	ext_cov_file = list()
	for line in coverage_file:
		line = line.strip('\n').split('\t')
		line.append(str(auto_mean_cov))
		line.append(str(total_mean_cov))
		line.append(str(gene_mean_cov[line[3].strip()]))
		if float(line[5].strip()) != 0:
			line.append(str(float(line[4].strip())/float(line[5].strip())))
		else:
			line.append(str(0))

		if float(line[6].strip()) != 0:
			line.append(str(float(line[4].strip())/float(line[6].strip())))
		else:
			line.append(str(0))
		
		if float(line[7].strip()) != 0:
			line.append(str(float(line[4].strip())/float(line[7].strip())))
		else:
			line.append(str(0))

		line[0] = line[0].strip()[3:]
		ext_cov_file.append(line)

	#Creating output path for saving
	output_path = StartWithMatchScore + '/' + get_sample_name(coverage_file_path) + '.aligned.only.normalized.coverage.txt'


	#Appending header at the beginning
	ext_cov_file.insert(0, ext_cov_file_header)

	#Storing extended coverage file
	with open(output_path, 'w') as file:
   		file.writelines('\t'.join(i) + '\n' for i in ext_cov_file)





	
	
