#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import numpy as np

from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import AdaBoostClassifier
def AB_AAC(train_name, train_label, test_name, test_label, AAC_Feature):
    second_train_label_feature = []
    for second_i in range(train_label.shape[0]):
        pair = []
        if train_label[second_i] == "1":  # continue
            continue
        if train_label[second_i] == "2":  # VF
            pair.append(1)
        if train_label[second_i] == "3":  # Negtive
            pair.append(0)
        for second_i3 in range(len(AAC_Feature)):
            if train_name[second_i] == AAC_Feature[second_i3][0].split("|")[1]:
                pair.extend(AAC_Feature[second_i3][1:])
                break
        second_train_label_feature.append(pair)
    second_test_label_feature = []
    for second_i in range(test_label.shape[0]):
        pair = []
        if test_label[second_i] == "1":
            pair.append(0)
        if test_label[second_i] == "2":
            pair.append(1)
        if test_label[second_i] == "3":
            pair.append(0)
        for second_i3 in range(len(AAC_Feature)):
            if test_name[second_i] == AAC_Feature[second_i3][0].split("|")[1]:
                pair.extend(AAC_Feature[second_i3][1:])
                break
        second_test_label_feature.append(pair)
    #########
    #########

    second_train_label_feature = np.array(second_train_label_feature, dtype=float)
    second_train_data = second_train_label_feature[:, 1:]
    second_train_label = second_train_label_feature[:, 0]
    second_train_label = second_train_label.astype(int)
    second_train_data = second_train_data.astype(float)
    second_test_label_feature = np.array(second_test_label_feature, dtype=float)
    second_test_data = second_test_label_feature[:, 1:]
    second_test_label = second_test_label_feature[:, 0]
    second_test_label = second_test_label.astype(int)
    second_test_data = second_test_data.astype(float)

    print('Train Data：', second_train_data.shape)
    print('Test Data：', second_test_data.shape)

    clf1 = AdaBoostClassifier(n_estimators=500)
    final_blend_train = []
    fianl_blend_test = []
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=10)
    for blend_train_index, blend_test_index in skf.split(second_train_data, second_train_label):
        blend_train_data, blend_train_label = second_train_data[blend_train_index], second_train_label[
            blend_train_index]
        blend_test_data, blend_test_label = second_train_data[blend_test_index], second_train_label[blend_test_index]
        clf1.fit(blend_train_data, blend_train_label)
        clf1_pred_proba = clf1.predict_proba(blend_test_data)[:, 1]
        for num_i in range(len(blend_test_label)):
            pair = []
            pair.append(blend_test_label[num_i])
            pair.append(clf1_pred_proba[num_i])
            final_blend_train.append(pair)
        ##############
        clf1_pred_proba2 = clf1.predict_proba(second_test_data)[:, 1]
        for num_j in range(len(second_test_label)):
            pair = []
            pair.append(second_test_label[num_j])
            pair.append(clf1_pred_proba2[num_j])
            fianl_blend_test.append(pair)
    fianl_blend_test_mean = []
    for num_k in range(int(len(fianl_blend_test) / 5)):
        pair = []
        pair.append(fianl_blend_test[num_k][0])
        for num_q in range(1, len(fianl_blend_test[num_k])):
            mean_num = (fianl_blend_test[num_k][num_q] + fianl_blend_test[num_k + len(second_test_label)][num_q] +
                        fianl_blend_test[num_k + len(second_test_label) * 2][num_q] +
                        fianl_blend_test[num_k + len(second_test_label) * 3][num_q] +
                        fianl_blend_test[num_k + len(second_test_label) * 4][num_q]) / 5
            pair.append(mean_num)
        fianl_blend_test_mean.append(pair)
    y_score0 = []
    for i in range(len(fianl_blend_test_mean)):
        if float(fianl_blend_test_mean[i][1]) >= 0.5:
            y_score0.append(1)
        else:
            y_score0.append(0)
    second_test_label = list(second_test_label)

    return final_blend_train, fianl_blend_test_mean