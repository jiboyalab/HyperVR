import os
import csv

from collections import defaultdict
import argparse
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(data)
    return


def GenerateBitscoreFeatureStep1(Remainingfilename):

    card_and_ardb_gene_name = open("Database_Gene_Name.txt", 'w')
    for i in open(Remainingfilename):
        i = i.strip().split(">")
        try:
            if i[1]:
                card_and_ardb_gene_name.writelines(i[1])
                card_and_ardb_gene_name.write("\n")
        except:
            continue
def GenerateBitscoreFeatureStep2(filename):
    final_pair = []
    with open("Database_Gene_Name.txt") as f:
        all_line = f.readlines()
        num = len(all_line)
    for i in open(filename):
        i = i.strip().split(">")
        pair = []
        try:
            if i[1]:
                pair.append(i[1])
                for k in range(num):
                    pair.append(0.0)
            final_pair.append(pair)
        except:
            continue
    StorFile(final_pair, "Feature0_Matrix.txt")
def GenerateBitscoreFeatureStep3(filename, Remainingfilename,Diamond_Path):

    os.system(Diamond_Path+ " makedb --in " + Remainingfilename + " --db Database_GENE.dmnd")
    os.system(Diamond_Path+
        " blastp --db Database_GENE.dmnd --query " + filename + " --out ArgVFNeg_diamond_Database.txt --more-sensitive")

    # with open("ArgVFNeg_diamond_Database.txt") as fr:
    #     all_lines = fr.readlines()
    # ArgVFNeg_diamond_Database2 = open("ArgVFNeg_diamond_Database2.txt", 'w')
    # for i in range(len(all_lines)):
    #     if all_lines[i].strip().split("\t")[0] != all_lines[i].strip().split("\t")[1]:
    #         ArgVFNeg_diamond_Database2.writelines(all_lines[i])

def GenerateBitscoreFeatureStep4():
    bit_score_feature_file = open("Bit_Score_Feature_Matrix.txt", 'w')

    with open("Feature0_Matrix.txt") as f:
        all_line = f.readlines()
        all_num = len(all_line)


    with open("Database_Gene_Name.txt") as f3:
        all_line3 = f3.readlines()
    Database_Gene_Name_dict=defaultdict(list)
    for i in range(len(all_line3)):
        seq=all_line3[i].strip()
        Database_Gene_Name_dict[seq].append(i)


    with open("ArgVFNeg_diamond_Database2.txt") as y:
        match_lines = y.readlines()
    ArgVFNeg_diamond_Database_dict = defaultdict(list)
    for i in range(len(match_lines)):
        seq = match_lines[i].strip().split("\t")
        ArgVFNeg_diamond_Database_dict[seq[0]].append(seq[1])
        ArgVFNeg_diamond_Database_dict[seq[0]].append(seq[-1])

    countnum = 0
    for i in open("Feature0_Matrix.txt"):
        print("countnum", countnum)
        countnum += 1
        i = i.strip().split(",")
        gene_name = i[0]
        list1=ArgVFNeg_diamond_Database_dict[gene_name]
        for j in range(0,len(list1),2):
            posi=Database_Gene_Name_dict[list1[j]][0]
            i[posi]=list1[j+1]


        for m in range(len(i) - 1):
            bit_score_feature_file.writelines(i[m])
            bit_score_feature_file.write(",")
        bit_score_feature_file.writelines(i[-1])
        if countnum < all_num:
            bit_score_feature_file.write("\n")


def GenerateBitscoreFeatureStep5():
    normal_bit_score_feature_file = open("Normalized_Bit_Score_Feature_Matrix.txt", 'w')
    with open("Bit_Score_Feature_Matrix.txt") as f:
        all_line = f.readlines()
        all_num = len(all_line)
    num = 0
    for i in open("Bit_Score_Feature_Matrix.txt"):
        num += 1
        print(num)
        i = i.strip().split(",")
        max = -1.0
        min = 10000.0
        normal_bit_score_feature_file.writelines(i[0])
        normal_bit_score_feature_file.write(",")
        for j in range(1, len(i)):
            # print j,i[j]
            now_num = float(i[j])
            if now_num < min:
                min = now_num
            if now_num > max:
                max = now_num
        if min == max == 0:
            for k in range(1, len(i) - 1):
                normal_bit_score_feature_file.writelines(str(i[k]))
                normal_bit_score_feature_file.write(",")
            normal_bit_score_feature_file.writelines(str(i[-1]))
        else:
            for k in range(1, len(i) - 1):
                normal_bit_score_feature_file.writelines(str((float(i[k]) - min) / (max - min)))
                normal_bit_score_feature_file.write(",")
            normal_bit_score_feature_file.writelines(str((float(i[-1]) - min) / (max - min)))
        if num < all_num:
            normal_bit_score_feature_file.write("\n")
def GenerateBitscoreFeature(filename, Remainingfilename,outputDir, Diamond_Path):

    GenerateBitscoreFeatureStep1(Remainingfilename)
    GenerateBitscoreFeatureStep2(filename)
    GenerateBitscoreFeatureStep3(filename, Remainingfilename,Diamond_Path)
    GenerateBitscoreFeatureStep4()
    GenerateBitscoreFeatureStep5()
    if os.path.exists(outputDir) == False:
        os.mkdir(outputDir)
    os.system("mv ArgVFNeg_diamond_Database.txt " + outputDir + "/Bitscore")
    os.system("mv ArgVFNeg_diamond_Database2.txt "+outputDir+"/Bitscore")
    os.system("mv Bit_Score_Feature_Matrix.txt " + outputDir + "/Bitscore")
    os.system("mv Database_GENE.dmnd " + outputDir + "/Bitscore")
    os.system("mv Database_Gene_Name.txt " + outputDir + "/Bitscore")
    os.system("mv Feature0_Matrix.txt " + outputDir + "/Bitscore")
    os.system("mv Normalized_Bit_Score_Feature_Matrix.txt " + outputDir + "/Bitscore")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage="it's usage tip.", description="generate bitscore feature")
    parser.add_argument("--file", required=True, help="protein sequence file in fasta format")
    parser.add_argument("--db_file", required=True, help="protein sequence file in fasta format")
    parser.add_argument("--diamond_path", help="the path of diamond program")
    parser.add_argument("--outdir", help="the path of out dir")
    args = parser.parse_args()

    filename=args.file
    Remainingfilename=args.db_file
    outputDir = args.outdir
    Diamond_Path=args.diamond_path

    GenerateBitscoreFeature(filename, Remainingfilename,outputDir, Diamond_Path)


