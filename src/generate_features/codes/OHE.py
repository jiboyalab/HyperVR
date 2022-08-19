#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import numpy as np

def OHE(fastas, **kw):
    # function to one-hot encode a Protein sequence string
    fixed_len=2000
    AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
    # AA = 'ARNDCQEGHILKMFPSTWYV'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        from numpy import argmax

        # define universe of possible input values
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        # define a mapping of chars to integers
        char_to_int = dict((c, i) for i, c in enumerate(alphabet))
        int_to_char = dict((i, c) for i, c in enumerate(alphabet))
        # integer encode input data
        integer_encoded = [char_to_int[char] for char in sequence]

        # one hot encode
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(len(alphabet))]
            letter[value] = 1
            onehot_encoded.append(letter)
        pair=[]
        if len(onehot_encoded) >= fixed_len:
            for i in range(fixed_len):
                pair.extend(onehot_encoded[i][:])
        else:
            for i in range(len(onehot_encoded)):
                pair.extend(onehot_encoded[i][:])
            for i in range(len(onehot_encoded),fixed_len):
                pair.extend([0]*20)
        code = code + pair
        encodings.append(code)

    return encodings


