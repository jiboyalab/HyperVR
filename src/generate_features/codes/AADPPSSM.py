#!/usr/bin/env python
#_*_coding:utf-8_*_
import numpy as np
import sys, os
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
import checkFasta
from AADPPSSMTools import *
def AADPPSSM(fastas, **kw):
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
		name, sequence = i[0], i[1]
		print(name)
		filename=i[0].replace('|','')
		length=len(sequence)
		code = [name]
		if os.path.exists(pssmDir+'/'+filename+'.pssm') == False:
			print('Error: pssm prfile for protein ' + name + ' does not exist.')#序列太短了
			temlist = ["0"] *420
			code=code+temlist
			encodings.append(code)
			continue
			# sys.exit(1)
		with open(pssmDir+'/'+filename+'.pssm') as f:
			records = f.readlines()[3: -6]
		proteinSeq = ''
		pssmMatrix = []
		for line in records:
			array = line.strip().split()
			pssmMatrix.append(array[1:42])
			proteinSeq = proteinSeq + array[1]

		pos = proteinSeq.find(sequence)
		if pos == -1:
			print('Warning: could not find the peptide in proteins.\n\n')
			temlist = ["0"] * 420
			code = code + temlist
			encodings.append(code)
			continue
		else:
			pair = []
			aadppssmMatrix=aadp_pssm(np.array(pssmMatrix))
			pair=pair+list(aadppssmMatrix[0])
			pair = [str(pair[i]) for i in range(len(pair))]
			code = code + pair
		encodings.append(code)


	return encodings