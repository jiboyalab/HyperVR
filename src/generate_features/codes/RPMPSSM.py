#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta

def RPMPSSM(fastas, **kw):
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
		# name, sequence = i[0], i[1]
		# length=len(sequence)
		# code = [name]
		if os.path.exists(pssmDir+'/'+filename+'.pssm') == False:
			print('Error: pssm prfile for protein ' + name + ' does not exist.')
			temlist = ["0"] * 400
			code = code + temlist
			encodings.append(code)
			continue
			# sys.exit(1)
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

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
		else:
			pair = []
			for aa in AA:
				temlist = [0] * 20
				for p in range(pos, pos + len(sequence)):
					aa2 = proteinSeq[p]
					aa2feature = pssmMatrix[p]
					if aa2 == aa:
						temlist = [temlist[i] + float(aa2feature[i]) for i in range(min(len(temlist), len(aa2feature)))]
				temlist = [(temlist[i]) / length for i in range(len(temlist))]
				pair = pair + temlist
			pair = [str(pair[i]) for i in range(len(pair))]
			code = code + pair
		encodings.append(code)


	return encodings