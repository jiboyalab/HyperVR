#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta

def SFPSSM(fastas, **kw):
	# if checkFasta.checkFasta(fastas) == False:
	# 	print('Error: for "PSSM" encoding, the input fasta sequences should be with equal length. \n\n')
	# 	return 0

	pssmDir = kw['path']
	if pssmDir == None:
		print('Error: please specify the directory of predicted protein disorder files by "--path" \n\n')
		return 0

	AA = 'ARNDCQEGHILKMFPSTWYV'

	encodings = []
	# header = ['#']
	# for p in range(1, len(fastas[0][1]) + 1):
	# 	for aa in AA:
	# 		header.append('Pos.'+str(p) + '.' + aa)
	# encodings.append(header)

	for i in fastas:
		name, sequence = i[0].split('|')[1], i[1]
		filename = i[0].replace('|', '')
		length = len(sequence)
		code = [name]

		if os.path.exists(pssmDir+'/'+filename+'.pssm') == False:
			print('Error: pssm prfile for protein ' + name + ' does not exist.')
			sys.exit(1)
		with open(pssmDir+'/'+filename+'.pssm') as f:
			records = f.readlines()[3: -6]
		proteinSeq = ''
		pssmMatrix = []
		for line in records:
			array = line.strip().split()
			pssmMatrix.append(array[2:22])
			proteinSeq = proteinSeq + array[1]

		for i in range(len(pssmMatrix)):
			for j in range(len(pssmMatrix[i])):
				if float(pssmMatrix[i][j])<0:
					pssmMatrix[i][j]="0"
				if float(pssmMatrix[i][j])>7:
					pssmMatrix[i][j] = "7"

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			pair=[]
			for aa in AA:
				num=0
				while num<20:
					sumnum = 0
					for p in range(pos, pos + len(sequence)):
						aa2=proteinSeq[p]
						aa2feature=pssmMatrix[p][num]
						if aa2==aa:
							sumnum+=float(aa2feature)
					pair.append(sumnum)
					num+=1
			pair=[str(pair[i]) for i in range(len(pair))]
			code = code+pair
		encodings.append(code)


	return encodings