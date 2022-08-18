#!/usr/bin/env python
#_*_coding:utf-8_*_

import sys, os,csv
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return
pPath = os.path.split(os.path.realpath(__file__))[0]
sys.path.append(pPath)
def MyKmer(MySequence, k=3):         # k是几mer，s是取前几个奇异值
    # 建立空KmerRow
    KmerRow = []
    counter = 0
    while counter < 20 ** k:
        KmerRow.append(0)
        counter = counter + 1
    # 向矩阵填1
    counter = 0
    while counter < len(MySequence) - k + 1:
        sequence = MySequence[counter:counter + k]
        # print(sequence)
        num = 0
        n = 0
        while n < k:
            if sequence[k - n - 1] == 'A':
                num = num + 0 * 20 ** n
            if sequence[k - n - 1] == 'R':
                num = num + 1 * 20 ** n
            if sequence[k - n - 1] == 'N':
                num = num + 2 * 20 ** n
            if sequence[k - n - 1] == 'D':
                num = num + 3 * 20 ** n
            if sequence[k - n - 1] == 'C':
                num = num + 4 * 20 ** n
            if sequence[k - n - 1] == 'Q':
                num = num + 5 * 20 ** n
            if sequence[k - n - 1] == 'E':
                num = num + 6 * 20 ** n
            if sequence[k - n - 1] == 'G':
                num = num + 7 * 20 ** n
            if sequence[k - n - 1] == 'H':
                num = num + 8 * 20 ** n
            if sequence[k - n - 1] == 'I':
                num = num + 9 * 20 ** n
            if sequence[k - n - 1] == 'L':
                num = num + 10 * 20 ** n
            if sequence[k - n - 1] == 'K':
                num = num + 11 * 20 ** n
            if sequence[k - n - 1] == 'M':
                num = num + 12 * 20 ** n
            if sequence[k - n - 1] == 'F':
                num = num + 13 * 20 ** n
            if sequence[k - n - 1] == 'P':
                num = num + 14 * 20 ** n
            if sequence[k - n - 1] == 'S':
                num = num + 15 * 20 ** n
            if sequence[k - n - 1] == 'T':
                num = num + 16 * 20 ** n
            if sequence[k - n - 1] == 'W':
                num = num + 17 * 20 ** n
            if sequence[k - n - 1] == 'Y':
                num = num + 18 * 20 ** n
            if sequence[k - n - 1] == 'V':
                num = num + 19 * 20 ** n

            n = n + 1
        # print(num)
        KmerRow[num] = KmerRow[num] + 1
        counter = counter + 1
    counter = 0
    while counter < len(KmerRow):
        KmerRow[counter] = KmerRow[counter] / (len(MySequence) - k + 1)
        counter = counter + 1
    return KmerRow

def kmer(fastas, **kw):
    encodings = []
    for i in fastas:
        pair=[]
        name, sequence = i[0].split('|')[1], i[1]
        feature=MyKmer(sequence)
        pair.append(name)
        pair.extend(feature)
        encodings.append(pair)
    StorFile(encodings, '/home/jby/PredVFandARG/UseOnlyKmer/KmerFeature.csv')
    return 0





