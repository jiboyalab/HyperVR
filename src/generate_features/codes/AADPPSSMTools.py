import numpy as np
def average(matrixSum, seqLen):
    # average the summary of rows
    matrix_array = np.array(matrixSum)
    matrix_array = np.divide(matrix_array, seqLen)
    matrix_array_shp = np.shape(matrix_array)
    matrix_average = [(np.reshape(matrix_array, (matrix_array_shp[0] * matrix_array_shp[1], )))]
    return matrix_average
def preHandleColumns(PSSM,STEP,PART,ID):
    '''
    if STEP=k, we calculate the relation betweem one residue and the kth residue afterward.
    '''
    '''
    if PART=0, we calculate the left part of PSSM.
    if PART=1, we calculate the right part of PSSM.
    '''
    '''
    if ID=0, we product the residue-pair.
    if ID=1, we minus the residue-pair.
    '''
    '''
    if KEY=1, we divide each element by the sum of elements in its column.
    if KEY=0, we don't perform the above process.
    '''
    if PART==0:
        #print "PART=",PART
        PSSM=PSSM[:,1:21]
    elif PART==1:
        #print "PART=",PART
        PSSM=PSSM[:, 21:]
    PSSM=PSSM.astype(float)
    matrix_final = [ [0.0] * 20 ] * 20
    matrix_final=np.array(matrix_final)
    seq_cn=np.shape(PSSM)[0]
    #print "seq_cn=",seq_cn

    if ID==0:
        #print "ID=",ID
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j]+=(PSSM[k][i]*PSSM[k+STEP][j])

    elif ID==1:
        #print "ID=",ID
        for i in range(20):
            for j in range(20):
                for k in range(seq_cn - STEP):
                    matrix_final[i][j] += ((PSSM[k][i]-PSSM[k+STEP][j]) * (PSSM[k][i]-PSSM[k+STEP][j])/4.0)
    #print matrix_final
    return matrix_final
def handleRows(PSSM, SWITCH, COUNT):
    '''
    if SWITCH=0, we filter no element.
    if SWITCH=1, we filter all the negative elements.
    if SWITCH=2, we filter all the negative and positive elements greater than expected.
    '''
    '''
    if COUNT=20, we generate a 20-dimension vector.
    if COUNT=400, we generate a 400-dimension vector.
    '''
    # 0-19 represents amino acid 'ARNDCQEGHILKMFPSTWYV'
    Amino_vec = "ARNDCQEGHILKMFPSTWYV"

    matrix_final = [ [0] * 20 ] * int(COUNT/20)
    # matrix_final=np.array(matrix_final)
    seq_cn = 0

    PSSM_shape=np.shape(PSSM)
    for i in range(PSSM_shape[0]):
        seq_cn += 1
        str_vec=PSSM[i]
        str_vec_positive=map(int, str_vec[1:21])
        # str_vec_positive=np.array(str_vec_positive)
        if SWITCH==1:
            str_vec_positive[str_vec_positive<0]=0
        elif SWITCH==2:
            str_vec_positive[str_vec_positive<0]=0
            str_vec_positive[str_vec_positive>7]=0
        #print "str_vec_positive="
        #print str_vec_positive
        if COUNT==20:
            k=0
            for u, i in zip(str_vec_positive, matrix_final[0]):
                matrix_final[0][k]=u+i
                k+=1

        elif COUNT==400:
            k = 0
            for u, i in zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])]):
                matrix_final[Amino_vec.index(str_vec[0])][k] = u + i
                k += 1
            # matrix_final[Amino_vec.index(str_vec[0])] = map(sum, zip(str_vec_positive, matrix_final[Amino_vec.index(str_vec[0])]))

        #print "matrix_final="
        #print matrix_final

    return matrix_final
def dpc_pssm(input_matrix):
    #print "start dpc_pssm function"
    PART = 0
    STEP = 1
    ID = 0
    KEY = 0
    matrix_final = preHandleColumns(input_matrix, STEP, PART, ID)
    seq_cn = float(np.shape(input_matrix)[0])
    dpc_pssm_vector = average(matrix_final, seq_cn-STEP)
    #print "end dpc_pssm function"
    return dpc_pssm_vector
def aac_pssm(input_matrix):
    #print "start aac_pssm function"
    SWITCH = 0
    COUNT = 20
    seq_cn=float(np.shape(input_matrix)[0])
    aac_pssm_matrix=handleRows(input_matrix,SWITCH,COUNT)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    aac_pssm_vector=average(aac_pssm_matrix,seq_cn)
    #print "end aac_pssm function"
    return aac_pssm_vector
def aadp_pssm(input_matrix):
    aac_pssm_matrix=aac_pssm(input_matrix)
    dpc_pssm_matrix=dpc_pssm(input_matrix)
    aac_pssm_matrix=np.array(aac_pssm_matrix)
    dpc_pssm_matrix=np.array(dpc_pssm_matrix)
    aadp_pssm_matrix=np.hstack((aac_pssm_matrix, dpc_pssm_matrix))
    #print np.shape(aadp_pssm_matrix)
    return aadp_pssm_matrix