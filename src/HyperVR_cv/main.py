#!/usr/bin/env python
#_*_coding:utf-8_*_
import argparse
import re,csv,math
from codes import *
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import ExtraTreesClassifier
import numpy as np
def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return
def ReadMyTsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName),delimiter = '\t')
    for row in csv_reader:  # 把每个rna疾病对加入OriginalData，注意表头
        SaveList.append(row)
    return
def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return
def MyConfusionMatrix(y_real,y_predict):
    from sklearn.metrics import confusion_matrix
    CM = confusion_matrix(y_real, y_predict)
    print(CM)
    CM = CM.tolist()
    TN = CM[0][0]
    FP = CM[0][1]
    FN = CM[1][0]
    TP = CM[1][1]
    print('TN:%d, FP:%d, FN:%d, TP:%d' % (TN, FP, FN, TP))
    Acc = (TN + TP) / (TN + TP + FN + FP)
    Sen = TP / (TP + FN)
    Spec = TN / (TN + FP)
    Prec = TP / (TP + FP)
    MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    Rec = TP/(TP+FN)
    F1=(2*Prec*Rec)  / (Prec+Rec)
    # 分母可能出现0，需要讨论待续
    # print('Acc:', round(Acc, 4))
    # print('Sen:', round(Sen, 4))
    # print('Spec:', round(Spec, 4))
    # print('Prec:', round(Prec, 4))
    # print('MCC:', round(MCC, 4))
    Result = []
    Result.append(round(Acc, 8))
    Result.append(round(Sen, 8))
    Result.append(round(Spec, 8))
    Result.append(round(Prec, 8))
    Result.append(round(MCC, 8))
    Result.append(round(Rec, 8))
    Result.append(round(F1, 8))
    return Result
def MyRealAndPredictionProb(Real,prediction):
    RealAndPredictionProb = []
    counter = 0
    while counter < len(prediction):
        pair = []
        pair.append(Real[counter])
        pair.append(prediction[counter][0])
        pair.append(prediction[counter][1])
        RealAndPredictionProb.append(pair)
        counter = counter + 1
    return RealAndPredictionProb
def MyRealAndPrediction(Real,prediction):
    RealAndPrediction = []
    counter = 0
    while counter < len(prediction):
        pair = []
        pair.append(Real[counter])
        pair.append(prediction[counter])
        RealAndPrediction.append(pair)
        counter = counter + 1
    return RealAndPrediction
if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.", description="Training and validating model by 5-fold cv")
    parser.add_argument("--feature_path", required=True, help="feature file path")
    parser.add_argument("--label_path", required=True, help="label file path")
    args = parser.parse_args()
    feature_path=args.feature_path
    label_path=args.label_path
    print("---please do not change the feature file name && Read features---\n")
    AAC_Feature = []
    ReadMyTsv(AAC_Feature, feature_path + "/AAC_encoding.tsv")
    AAC_Feature = AAC_Feature[1:][:] # Delete table header
    DDE_Feature = []
    ReadMyTsv(DDE_Feature, feature_path + "/DDE_encoding.tsv")
    DDE_Feature = DDE_Feature[1:][:]  # Delete table header
    DPC_Feature = []
    ReadMyTsv(DPC_Feature, feature_path + "/DPC_encoding.tsv")
    DPC_Feature = DPC_Feature[1:][:]  # Delete table header
    PAAC_Feature = []
    ReadMyTsv(PAAC_Feature,feature_path+"/PAAC_encoding.tsv")
    PAAC_Feature = PAAC_Feature[1:][:] # Delete table header
    QSO_Feature = []
    ReadMyTsv(PAAC_Feature, feature_path + "/QSOrder_encoding.tsv")
    QSO_Feature = QSO_Feature[1:][:] # Delete table header
    OneHot_Feature = []
    ReadMyTsv(OneHot_Feature, feature_path + "/OHE_encoding.tsv")
    OneHot_Feature = OneHot_Feature[1:][:] # Delete table header
    Bit_Score_Feature = []
    ReadMyCsv(Bit_Score_Feature, feature_path + "/bitscore/Normalized_Bit_Score_Feature_Matrix.txt")#
    PSSM_Com_Feature = []
    ReadMyTsv(PSSM_Com_Feature, feature_path + "/PSSMC_encoding.tsv")
    AADPPSSM_Feature = []
    ReadMyTsv(AADPPSSM_Feature, feature_path + "/AADPPSSM_encoding.tsv")
    RPMPSSM_Feature = []
    ReadMyTsv(RPMPSSM_Feature, feature_path + "/RPMPSSM_encoding.tsv")
    Label=[]
    ReadMyCsv(Label, label_path + "/Label_ARG+VF+NS.csv.csv")#
    Label = np.array(Label)
    name = Label[:, 0]
    label = Label[:, 1]
    print("---start cross validation analysis---\n")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=10)
    k=0
    for train_index, test_index in skf.split(name, label):
        train_name, train_label = name[train_index], label[train_index]
        test_name, test_label = name[test_index], label[test_index]
        print("--HyperVR-ARG---\n")
        final_blend_train = []
        fianl_blend_test_mean = []
        print("---deep learning && bitscore features---\n")
        final_blend_train1, fianl_blend_test_mean1 = DL_Bitscore.DL_Bitscore(train_name, train_label, test_name,test_label,Bit_Score_Feature)
        print("---extratree classifier && PSSMC features---\n")
        final_blend_train2, fianl_blend_test_mean2 = ET_Pssmc.ET_Pssmc(train_name, train_label, test_name,test_label, PSSM_Com_Feature)
        print("---extratree classifier && AADPPSSM features---\n")
        final_blend_train3, fianl_blend_test_mean3 = ET_AADPPSSM.ET_AADPPSSM(train_name, train_label, test_name, test_label,AADPPSSM_Feature)
        print("---extratree classifier && RPMPPSSM features---\n")
        final_blend_train4, fianl_blend_test_mean4 = ET_RPMPSSM.ET_RPMPSSM(train_name, train_label, test_name, test_label, RPMPSSM_Feature)
        print("---Stacking model && All features---\n")
        for i in range(len(final_blend_train1)):
            pair = []
            pair.append(final_blend_train1[i][0])
            pair.append(final_blend_train1[i][1])
            pair.append(final_blend_train2[i][1])
            pair.append(final_blend_train3[i][1])
            pair.append(final_blend_train4[i][1])
            final_blend_train.append(pair)
        for i in range(len(fianl_blend_test_mean1)):
            pair = []
            pair.append(fianl_blend_test_mean1[i][0])
            pair.append(fianl_blend_test_mean1[i][1])
            pair.append(fianl_blend_test_mean2[i][1])
            pair.append(fianl_blend_test_mean3[i][1])
            pair.append(fianl_blend_test_mean4[i][1])
            fianl_blend_test_mean.append(pair)
        final_blend_train = np.array(final_blend_train, dtype=float)
        fianl_blend_test_mean = np.array(fianl_blend_test_mean, dtype=float)
        final_blend_train_data, final_blend_train_label = final_blend_train[:, 1:], final_blend_train[:, 0]
        final_blend_test_data, final_blend_test_label = fianl_blend_test_mean[:, 1:], fianl_blend_test_mean[:, 0]
        clfStack = ExtraTreesClassifier(n_estimators=500)
        clfStack.fit(final_blend_train_data, final_blend_train_label)
        y_score0 = clfStack.predict(final_blend_test_data)
        y_score1 = clfStack.predict_proba(final_blend_test_data)
        RealAndPrediction = MyRealAndPrediction(final_blend_test_label, y_score0) #Real and predicted labels for the first step ARG prediction
        RealAndPredictionProb = MyRealAndPredictionProb(final_blend_test_label, y_score1) #Real and predicted probabilities for the first step ARG prediction
        StorFile(RealAndPrediction,str(k)+"RealAndPredictionForARG.csv") # k fold
        StorFile(RealAndPredictionProb, str(k)+"RealAndPredictionProbForARG.csv")# k fold
        clfStackResult = MyConfusionMatrix(final_blend_test_label, y_score0)
        print("Stacking model (HyperVR-ARG) Test accuracy: precision: recall: f1:", clfStackResult[0], clfStackResult[3],
              clfStackResult[5], clfStackResult[6])
        print("--HyperVR-VF---\n")
        test_name2 = []
        test_label2 = []

        for i in range(len(RealAndPrediction)):# Predicting the remaining genes that are not ARGs
            if RealAndPrediction[i][1] == 0.0:
                test_name2.append(name[test_index[i]])
                test_label2.append(label[test_index[i]])
        test_name = np.array(test_name2)
        test_label = np.array(test_label2)
        final_blend_train = []
        fianl_blend_test_mean = []
        print("---random forest && AAC features---\n")
        final_blend_train1, fianl_blend_test_mean1 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, AAC_Feature)
        print("---xgboost && AAC features---\n")
        final_blend_train2, fianl_blend_test_mean2 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, AAC_Feature)
        print("---extratree classifier && AAC features---\n")
        final_blend_train3, fianl_blend_test_mean3 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label, AAC_Feature)
        print("---gradientboosting && AAC features---\n")
        final_blend_train4, fianl_blend_test_mean4 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label, AAC_Feature)
        print("---adaboost && AAC features---\n")
        final_blend_train5, fianl_blend_test_mean5 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, AAC_Feature)

        print("---random forest && DPC features---\n")
        final_blend_train6, fianl_blend_test_mean6 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, DPC_Feature)
        print("---xgboost && DPC features---\n")
        final_blend_train7, fianl_blend_test_mean7 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, DPC_Feature)
        print("---extratree classifier && DPC features---\n")
        final_blend_train8, fianl_blend_test_mean8 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label, DPC_Feature)
        print("---gradientboosting && DPC features---\n")
        final_blend_train9, fianl_blend_test_mean9 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label,  DPC_Feature)
        print("---adaboost && DPC features---\n")
        final_blend_train10, fianl_blend_test_mean10 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, DPC_Feature)

        print("---random forest && DDE features---\n")
        final_blend_train11, fianl_blend_test_mean11 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, DDE_Feature)
        print("---xgboost && DDE features---\n")
        final_blend_train12, fianl_blend_test_mean12 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name,  test_label, DDE_Feature)
        print("---extratree classifier && DDE features---\n")
        final_blend_train13, fianl_blend_test_mean13 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label,  DDE_Feature)
        print("---gradientboosting && DDE features---\n")
        final_blend_train14, fianl_blend_test_mean14 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label,  DDE_Feature)
        print("---adaboost && DDE features---\n")
        final_blend_train15, fianl_blend_test_mean15 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, DDE_Feature)

        print("---random forest && PAAC features---\n")
        final_blend_train16, fianl_blend_test_mean16 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, PAAC_Feature)
        print("---xgboost && PAAC features---\n")
        final_blend_train17, fianl_blend_test_mean17 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, PAAC_Feature)
        print("---extratree classifier && PAAC features---\n")
        final_blend_train18, fianl_blend_test_mean18 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label,  PAAC_Feature)
        print("---gradientboosting && PAAC features---\n")
        final_blend_train19, fianl_blend_test_mean19 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label, PAAC_Feature)
        print("---adaboost && PAAC features---\n")
        final_blend_train20, fianl_blend_test_mean20 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, PAAC_Feature)

        print("---random forest && QSO features---\n")
        final_blend_train21, fianl_blend_test_mean21 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, QSO_Feature)
        print("---xgboost && QSO features---\n")
        final_blend_train22, fianl_blend_test_mean22 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, QSO_Feature)
        print("---extratree classifier && QSO features---\n")
        final_blend_train23, fianl_blend_test_mean23 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label,  QSO_Feature)
        print("---gradientboosting && QSO features---\n")
        final_blend_train24, fianl_blend_test_mean24 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label, QSO_Feature)
        print("---adaboost && QSO features---\n")
        final_blend_train25, fianl_blend_test_mean25 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, QSO_Feature)

        print("---random forest && PSSMC features---\n")
        final_blend_train26, fianl_blend_test_mean26 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, PSSM_Com_Feature)
        print("---xgboost && PSSMC features---\n")
        final_blend_train27, fianl_blend_test_mean27 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, PSSM_Com_Feature)
        print("---extratree classifier && PSSMC features---\n")
        final_blend_train28, fianl_blend_test_mean28 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label,  PSSM_Com_Feature)
        print("---gradientboosting && PSSMC features---\n")
        final_blend_train29, fianl_blend_test_mean29 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label, PSSM_Com_Feature)
        print("---adaboost && PSSMC features---\n")
        final_blend_train30, fianl_blend_test_mean30 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, PSSM_Com_Feature)

        print("---random forest && AADPPSSM features---\n")
        final_blend_train31, fianl_blend_test_mean31 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label, AADPPSSM_Feature)
        print("---xgboost && AADPPSSM features---\n")
        final_blend_train32, fianl_blend_test_mean32 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, AADPPSSM_Feature)
        print("---extratree classifier && AADPPSSM features---\n")
        final_blend_train33, fianl_blend_test_mean33 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label,AADPPSSM_Feature)
        print("---gradientboosting && AADPPSSM features---\n")
        final_blend_train34, fianl_blend_test_mean34 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label, AADPPSSM_Feature)
        print("---adaboost && AADPPSSM features---\n")
        final_blend_train35, fianl_blend_test_mean35 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, AADPPSSM_Feature)

        print("---random forest && RPMPSSM features---\n")
        final_blend_train36, fianl_blend_test_mean36 = RF_AAC.RF_AAC(train_name, train_label, test_name, test_label,RPMPSSM_Feature)
        print("---xgboost && RPMPSSM features---\n")
        final_blend_train37, fianl_blend_test_mean37 = Xgboost_AAC.Xgboost_AAC(train_name, train_label, test_name, test_label, RPMPSSM_Feature)
        print("---extratree classifier && RPMPSSM features---\n")
        final_blend_train38, fianl_blend_test_mean38 = ET_AAC.ET_AAC(train_name, train_label, test_name, test_label, RPMPSSM_Feature)
        print("---gradientboosting && RPMPSSM features---\n")
        final_blend_train39, fianl_blend_test_mean39 = GB_AAC.GB_AAC(train_name, train_label, test_name, test_label,RPMPSSM_Feature)
        print("---adaboost && RPMPSSM features---\n")
        final_blend_train40, fianl_blend_test_mean40 = AB_AAC.AB_AAC(train_name, train_label, test_name, test_label, RPMPSSM_Feature)

        print("---deep learning && one hot features---\n")
        final_blend_train41, fianl_blend_test_mean41 = DL_OHE.DL_OHE(train_name, train_label, test_name,test_label, OneHot_Feature)

        print("---Stacking model && All features---\n")
        for i in range(len(final_blend_train1)):
            pair = []
            pair.append(final_blend_train1[i][0])
            pair.append(final_blend_train1[i][1])
            pair.append(final_blend_train2[i][1])
            pair.append(final_blend_train3[i][1])
            pair.append(final_blend_train4[i][1])
            pair.append(final_blend_train5[i][1])
            pair.append(final_blend_train6[i][1])
            pair.append(final_blend_train7[i][1])
            pair.append(final_blend_train8[i][1])
            pair.append(final_blend_train9[i][1])
            pair.append(final_blend_train10[i][1])
            pair.append(final_blend_train11[i][1])
            pair.append(final_blend_train12[i][1])
            pair.append(final_blend_train13[i][1])
            pair.append(final_blend_train14[i][1])
            pair.append(final_blend_train15[i][1])
            pair.append(final_blend_train16[i][1])
            pair.append(final_blend_train17[i][1])
            pair.append(final_blend_train18[i][1])
            pair.append(final_blend_train19[i][1])
            pair.append(final_blend_train20[i][1])
            pair.append(final_blend_train21[i][1])
            pair.append(final_blend_train22[i][1])
            pair.append(final_blend_train23[i][1])
            pair.append(final_blend_train24[i][1])
            pair.append(final_blend_train25[i][1])
            pair.append(final_blend_train26[i][1])
            pair.append(final_blend_train27[i][1])
            pair.append(final_blend_train28[i][1])
            pair.append(final_blend_train29[i][1])
            pair.append(final_blend_train30[i][1])
            pair.append(final_blend_train31[i][1])
            pair.append(final_blend_train32[i][1])
            pair.append(final_blend_train33[i][1])
            pair.append(final_blend_train34[i][1])
            pair.append(final_blend_train35[i][1])
            pair.append(final_blend_train36[i][1])
            pair.append(final_blend_train37[i][1])
            pair.append(final_blend_train38[i][1])
            pair.append(final_blend_train39[i][1])
            pair.append(final_blend_train40[i][1])
            pair.append(final_blend_train41[i][1])
            final_blend_train.append(pair)
        for i in range(len(fianl_blend_test_mean1)):
            pair = []
            pair.append(fianl_blend_test_mean1[i][0])
            pair.append(fianl_blend_test_mean1[i][1])
            pair.append(fianl_blend_test_mean2[i][1])
            pair.append(fianl_blend_test_mean3[i][1])
            pair.append(fianl_blend_test_mean4[i][1])
            pair.append(fianl_blend_test_mean5[i][1])
            pair.append(fianl_blend_test_mean6[i][1])
            pair.append(fianl_blend_test_mean7[i][1])
            pair.append(fianl_blend_test_mean8[i][1])
            pair.append(fianl_blend_test_mean9[i][1])
            pair.append(fianl_blend_test_mean10[i][1])
            pair.append(fianl_blend_test_mean11[i][1])
            pair.append(fianl_blend_test_mean12[i][1])
            pair.append(fianl_blend_test_mean13[i][1])
            pair.append(fianl_blend_test_mean14[i][1])
            pair.append(fianl_blend_test_mean15[i][1])
            pair.append(fianl_blend_test_mean16[i][1])
            pair.append(fianl_blend_test_mean17[i][1])
            pair.append(fianl_blend_test_mean18[i][1])
            pair.append(fianl_blend_test_mean19[i][1])
            pair.append(fianl_blend_test_mean20[i][1])
            pair.append(fianl_blend_test_mean21[i][1])
            pair.append(fianl_blend_test_mean22[i][1])
            pair.append(fianl_blend_test_mean23[i][1])
            pair.append(fianl_blend_test_mean24[i][1])
            pair.append(fianl_blend_test_mean25[i][1])
            pair.append(fianl_blend_test_mean26[i][1])
            pair.append(fianl_blend_test_mean27[i][1])
            pair.append(fianl_blend_test_mean28[i][1])
            pair.append(fianl_blend_test_mean29[i][1])
            pair.append(fianl_blend_test_mean30[i][1])
            pair.append(fianl_blend_test_mean31[i][1])
            pair.append(fianl_blend_test_mean32[i][1])
            pair.append(fianl_blend_test_mean33[i][1])
            pair.append(fianl_blend_test_mean34[i][1])
            pair.append(fianl_blend_test_mean35[i][1])
            pair.append(fianl_blend_test_mean36[i][1])
            pair.append(fianl_blend_test_mean37[i][1])
            pair.append(fianl_blend_test_mean38[i][1])
            pair.append(fianl_blend_test_mean39[i][1])
            pair.append(fianl_blend_test_mean40[i][1])
            pair.append(fianl_blend_test_mean41[i][1])
            fianl_blend_test_mean.append(pair)
        final_blend_train = np.array(final_blend_train, dtype=float)
        fianl_blend_test_mean = np.array(fianl_blend_test_mean, dtype=float)
        final_blend_train_data, final_blend_train_label = final_blend_train[:, 1:], final_blend_train[:, 0]
        final_blend_test_data, final_blend_test_label = fianl_blend_test_mean[:, 1:], fianl_blend_test_mean[:, 0]
        clfStack = ExtraTreesClassifier(n_estimators=500)
        clfStack.fit(final_blend_train_data, final_blend_train_label)
        y_score0 = clfStack.predict(final_blend_test_data)  #
        y_score1 = clfStack.predict_proba(final_blend_test_data)  #
        clfStackResult = MyConfusionMatrix(final_blend_test_label, y_score0)
        RealAndPrediction = MyRealAndPrediction(final_blend_test_label, y_score0)#Real and predicted labels for the Second step VF prediction
        RealAndPredictionProb = MyRealAndPredictionProb(final_blend_test_label, y_score1)#Real and predicted probabilities for the Second step VF prediction
        StorFile(RealAndPrediction, str(k) + "RealAndPredictionForVF.csv")  # k fold
        StorFile(RealAndPredictionProb, str(k) + "RealAndPredictionProbForVF.csv")  # k fold
        print("Stacking model (HyperVR-VF) Test accuracy: precision: recall: f1:", clfStackResult[0], clfStackResult[3],
              clfStackResult[5], clfStackResult[6])
        k+=1