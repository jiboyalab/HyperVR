#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os, re
pPath = re.sub(r'scripts$', '', os.path.split(os.path.realpath(__file__))[0])
sys.path.append(pPath)

import argparse

#!/usr/bin/env python
#_*_coding:utf-8_*_

import re, os, sys

def readFasta(file):
	if os.path.exists(file) == False:
		print('Error: "' + file + '" does not exist.')
		sys.exit(1)

	with open(file) as f:
		records = f.read()

	if re.search('>', records) == None:
		print('The input file seems not in fasta format.')
		sys.exit(1)

	records = records.split('>')[1:]
	myFasta = []
	for fasta in records:
		array = fasta.split('\n')
		name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV]', '', ''.join(array[1:]).upper())
		myFasta.append([name, sequence])
	return myFasta
def generatePSSMProfile(fastas, outDir, blastpgp, db):
	"""
	Generate PSSM file by using the blastpgp program in NCBI blast-2.2.18 package.

	Parameters
	----------
	file : file
		the file, which include the protein sequences in fasta format.

	blastpgp: string
		the path of blastpgp program.
		
	db: string 
		the uniref50 data base, which is formated by the 'formatdb' program in blast package.

	Returns
	-------
	a string: 
		A directory name, which include the predicted protein PSSM information.
	"""


	if os.path.exists(outDir) == False:
		os.mkdir(outDir)

	for i in fastas:
		name, sequence = re.sub('\|', '', i[0]), i[1]
		with open(name + '.txt', 'w') as f:
			f.write('>'+name+'\n'+sequence + '\n')
		myCmd = blastpgp + ' -query ' + name + '.txt' + ' -db ' + db + ' -num_threads 280 -num_iterations 3 -out_ascii_pssm ' + outDir + '/' + name +'.pssm'
		print('Doing psiblast for protein: ' + name)
		os.system(myCmd)
		os.remove(name + '.txt')
	return outDir


if __name__ == '__main__':
	"""
	Default path or specified by parameter
	"""
	dbName = '/home/jby/Blast/uniref50_blast'
	ncbidir = '/home/jby/Blast/ncbi-blast-2.12.0+/bin'
	outdirname = "/home/jby/Blast/outDir"
	parser = argparse.ArgumentParser(usage="it's usage tip.", description="generate PSSM profile")
	parser.add_argument("--file", required=True, help="protein sequence file in fasta format")
	parser.add_argument("--blastpgp", help="the path of NCBI psiblast program")
	parser.add_argument("--db", help="the path of unief50 database")
	parser.add_argument("--outdir", help="the path of out dir")
	args = parser.parse_args()

	blastpgp = args.blastpgp if args.blastpgp != None else ncbidir + '/psiblast'
	db = args.db if args.db != None else dbName
	outDir = args.outdir if args.db != None else outdirname
	fastas = readFasta(args.file)
	outputDir = generatePSSMProfile(fastas, outDir, blastpgp, db)
	print('The PSSM profiles are stored in directory: ' + outputDir)








