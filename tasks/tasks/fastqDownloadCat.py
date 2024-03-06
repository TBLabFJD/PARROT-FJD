#!/usr/bin/python

#########################################################################################
### Task:  Download samples from BASESPACE and/or concatenate FASTQ files for mapping ###
#########################################################################################

import sys
import os
from glob import glob
import subprocess
import argparse
import datetime
import shutil
import ConfigParser



def joinFastq(basespaceF, fastqFolder, inputF, dnaid, user, bscp):

	
	if basespaceF=="True":
		dnaid2 = dnaid
		#dnaid=dnaid[0:7]

		sampleFolder = fastqFolder+'/'+dnaid2
		os.mkdir(sampleFolder) # check if existing and delete content.
		sys.stdout.write("Downloading fastq files"+"\n")
		
		if user=="False":
			
			config = "."

		else:

			config = user


		cmd = bscp + " -q -s //"+ config +"/Projects/"+inputF+"/Samples/"+dnaid2+" "+sampleFolder

		#sys.stderr.write(cmd)
		proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		stdout, stderr = proc.communicate()
		exitcode = proc.returncode
		
		if exitcode != 0:
			sys.stderr.write("ERROR: '%s'\n" %(stderr))

			if "More than one matching identifier found for" in stderr:
				id = stderr.split("\n")[1].split('(', 1)[1].split(')')[0]
				sys.stderr.write("WARNING: Downloading sample '%s'\n" %(stderr.split("\n")[1].strip()))
				cmd= bscp + " -q -s //"+ config +"/Projects/"+inputF+"/Samples/"+id+" "+sampleFolder
				exitcode = subprocess.call(cmd, shell=True)
				if exitcode != 0:
					sys.stderr.write("ERROR: 'Something was wrong when downloading sample %s from basespace'.\n" %(stderr.split("\n")[1]))
			else:
				sys.exit(1)


		foward = sorted(glob(sampleFolder+ '/*_R1*.f*q.gz'))
		reverse = sorted(glob(sampleFolder+ '/*_R2*.f*q.gz'))

	else:
		foward = sorted(glob(inputF+ '/' + dnaid + '*_R1*.f*q.gz'))
		reverse = sorted(glob(inputF+ '/' + dnaid + '*_R2*.f*q.gz'))
	

	if len(foward)==len(reverse) and len(foward)!= 0:
		
		### FASTQ CONCATENATION	

		foward_file = fastqFolder + dnaid + '_R1.fastq.gz'
		reverse_file = fastqFolder + dnaid + '_R2.fastq.gz'

		sys.stdout.write(foward_file + " written"+"\n")
		sys.stdout.write(reverse_file +" written"+"\n")

		#sys.stdout.write('Joining forward fastq files of sample ' + dnaid +"\n")
		mycmd = 'cat %s > %s' %(' '.join(foward), foward_file)
		subprocess.call(mycmd, shell=True)


		#sys.stdout.write('Joining reverse fastq files of sample ' + dnaid +"\n")
		mycmd2 = 'cat %s > %s' %(' '.join(reverse), reverse_file)
		subprocess.call(mycmd2, shell=True)
		
	else:
		
		sys.stderr.write("ERROR: Not fastq files found for sample '%s' or different number of reverse and foward fastq files.\n" %(dnaid))
		sys.exit()


	if basespaceF=="True":
		shutil.rmtree(sampleFolder)

	# if basespaceF=="True":
	# 	shutil.rmtree(sampleFolder)


if __name__ == "__main__":

	# arguments
	# 1. basespace downloading?
	# 2. folder where store fastqs
	# 3. project to download or folder with fastq to concatenate
	# 4. sample name
	# 5. user

	configFilePath = os.path.dirname(os.path.realpath(__file__))+"/../pipeline.config"
	configFile = ConfigParser.ConfigParser()
	configFile.read(configFilePath)
	bscp=configFile.get("configFilePipeline","baseSpacecp_bin").strip('"').strip('\'')

	joinFastq(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], bscp)



