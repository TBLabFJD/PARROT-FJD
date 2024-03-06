#!/usr/bin/env python

import os
import argparse
import sys


def filtergenes(myargs):

	geneList = []

	with open(myargs.geneFilt, "r") as genes:

		for gene in genes:
			geneList.append(gene.strip().split("\t")[0])

	output = open(myargs.output, "w")

	with open(myargs.input, "r") as variants:

		first_line = variants.readline()
		output.write(first_line)


		for variant in variants:

			if variant.strip().split("\t")[11] in geneList:

				output.write(variant)


	output.close()






def main():

	#arguments
	parser = argparse.ArgumentParser(description="Fileting SNVs by gene list")
	parser.add_argument('-i', '--input', help='\t\tInput PVM file', required=True)
	parser.add_argument('-o', '--output', dest="output",help='\t\tOutput file.', required=True)
	parser.add_argument('-f', '--geneFilt', help='\t\tFinal filter of genes', required=True)

	
	args = parser.parse_args()


	# check if files and dirs exist
	
	if not os.path.isfile(args.input): 
		sys.stderr.write("ERROR: PVM file with called SNVs ('%s') does not exist\n" %(args.input))
		sys.exit()


	if not os.path.isfile(args.geneFilt): 
		sys.stderr.write("ERROR: File listing genes of interest ('%s') does not exist\n" %(args.geneFilt))
		sys.exit()


	dir_path = os.path.dirname(args.output)
	if not os.path.isdir(dir_path): 
		sys.stderr.write("ERROR: Output path ('%s') does not exist\n" %(args.output))
		sys.exit()



	filtergenes(args)



if __name__ == "__main__":
	main()







