#!/usr/bin/env python
#_*_coding:utf-8_*_

import re
import numpy as np
from sklearn.model_selection import StratifiedKFold
import tensorflow as tf
from tensorflow import keras
from tensorflow.python.keras.callbacks import EarlyStopping
def First_Model_DNN(blend_train_data,blend_train_label,blend_test_data,blend_test_label,second_test_data):
	input_dim = blend_train_data.shape[1]
	output_dim = len(set(blend_train_label))
	inputs = keras.Input(shape=(input_dim,))
	x = keras.layers.Dense(2 ** 12, activation="relu")(inputs)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	x = keras.layers.Dense(2 ** 10, activation="relu")(x)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	x = keras.layers.Dense(2 ** 8, activation="relu")(x)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	x = keras.layers.Dense(2 ** 6, activation="relu")(x)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	x = keras.layers.Dense(2 ** 4, activation="relu")(x)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	x = keras.layers.Dense(2 ** 2, activation="relu")(x)
	x = keras.layers.Dropout(0.05, noise_shape=None, seed=None)(x)
	outputs = keras.layers.Dense(1, activation="sigmoid")(x)
	model = keras.Model(inputs, outputs)
	model.compile(
		optimizer=keras.optimizers.SGD(learning_rate=0.05, momentum=0.9),
		loss=keras.losses.BinaryCrossentropy(),
		metrics=[keras.metrics.BinaryAccuracy(name="acc")],
	)
	batch_size = blend_train_data.shape[0]
	train_dataset = tf.data.Dataset.from_tensor_slices((blend_train_data, blend_train_label)).batch(batch_size)
	test_dataset = tf.data.Dataset.from_tensor_slices((blend_test_data, blend_test_label)).batch(batch_size)
	test_dataset2 = tf.data.Dataset.from_tensor_slices((second_test_data)).batch(batch_size)
	callbacks = [
		EarlyStopping('val_acc', patience=100),
		keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.01,
										  patience=100, min_lr=0.01)
	]
	history = model.fit(train_dataset, epochs=500, validation_data=test_dataset, callbacks=callbacks)

	# Name=str(k)+str(k2)+"First_model_DNN+Bitscore_Feature.txt"
	# np.save(Name, history.history)
	# Test the model on all available devices.
	prediction1 = model.predict(test_dataset)
	prediction2 = model.predict(test_dataset2)
	return prediction1, prediction2
def DL_Bitscore(train_name, train_label, test_name,test_label,Bit_Score_Feature):
	second_train_label_feature = []
	for second_i in range(train_label.shape[0]):
		pair = []
		if train_label[second_i] == "1":  # VF && continue
			pair.append(1)
		else:
			pair.append(0)
		for second_i3 in range(len(Bit_Score_Feature)):
			if train_name[second_i] == Bit_Score_Feature[second_i3][0]:
				pair.extend(Bit_Score_Feature[second_i3][1:])
				break
		second_train_label_feature.append(pair)
	second_test_label_feature = []
	for second_i in range(test_label.shape[0]):
		pair = []
		if test_label[second_i] == "1":
			pair.append(1)
		else:
			pair.append(0)
		for second_i3 in range(len(Bit_Score_Feature)):
			if test_name[second_i] == Bit_Score_Feature[second_i3][0]:
				pair.extend(Bit_Score_Feature[second_i3][1:])
				break
		second_test_label_feature.append(pair)

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
	print('Train data：', second_train_data.shape)
	print('Test Data：', second_test_data.shape)
	final_blend_train = []
	fianl_blend_test = []
	k2 = 0
	skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=10)
	for blend_train_index, blend_test_index in skf.split(second_train_data, second_train_label):
		blend_train_data, blend_train_label = second_train_data[blend_train_index], second_train_label[
			blend_train_index]
		blend_test_data, blend_test_label = second_train_data[blend_test_index], second_train_label[
			blend_test_index]
		clf1_pred_proba, clf1_pred_proba2 = First_Model_DNN(blend_train_data, blend_train_label, blend_test_data,
															blend_test_label, second_test_data)

		for num_i in range(len(blend_test_label)):
			pair = []
			pair.append(blend_test_label[num_i])
			pair.append(clf1_pred_proba[num_i])
			final_blend_train.append(pair)
		##############

		for num_j in range(len(second_test_label)):
			pair = []
			pair.append(second_test_label[num_j])
			pair.append(clf1_pred_proba2[num_j])
			fianl_blend_test.append(pair)
		k2 += 1
	fianl_blend_test_mean = []
	for num_k in range(int(len(fianl_blend_test) / 5)):  # Reconstruct the test data and take the average of 5 times
		pair = []
		pair.append(fianl_blend_test[num_k][0])
		for num_q in range(1, len(fianl_blend_test[num_k])):  # one column is label
			mean_num = (fianl_blend_test[num_k][num_q] + fianl_blend_test[num_k + len(second_test_label)][num_q] +
						fianl_blend_test[num_k + len(second_test_label) * 2][num_q] +
						fianl_blend_test[num_k + len(second_test_label) * 3][num_q] +
						fianl_blend_test[num_k + len(second_test_label) * 4][num_q]) / 5
			pair.append(mean_num[0])
		fianl_blend_test_mean.append(pair)
	y_score0 = []
	for i in range(len(fianl_blend_test_mean)):
		if float(fianl_blend_test_mean[i][1]) >= 0.5:
			y_score0.append(1)
		else:
			y_score0.append(0)
	second_test_label = list(second_test_label)


	return final_blend_train, fianl_blend_test_mean