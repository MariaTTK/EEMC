def evaluation(new_clusters,COMPLEXES_):
    #global COMPLEXES_
    import numpy as np
    import math

    np.seterr(divide='ignore', invalid='ignore')
    
    T=np.zeros((len(COMPLEXES_),len(new_clusters)))
    for i in range(len(COMPLEXES_)):
        if (not np.nonzero(COMPLEXES_[i]))==False:
            for j in range(len(new_clusters)):
               T[i,j]=len(list(set(COMPLEXES_[i]).intersection(new_clusters[j])))
        else:
               T[i,:]=0
               
        
    complex_sizes=np.zeros((len(COMPLEXES_),1))
    for i in range(len(COMPLEXES_)):
        complex_sizes[i]=len(COMPLEXES_[i])

    Sn=np.zeros((len(COMPLEXES_),len(new_clusters)))
    for j in range(len(new_clusters)):
        Sn[:,j]=T[:,j]/complex_sizes[j]
               
    Sn_co=Sn.max(1)    #complex wise sensitivity
                        #a.max(1),maximum element of each row of matrix a 
    for i in range(len(Sn_co)):
        if (math.isnan(Sn_co.tolist()[i])==True) or (math.isinf(Sn_co.tolist()[i])==True): #if itme=Nan then item=0
            Sn_co[i]=0

    #Sensitivity measures the proportion of actual positives which are correctly identified as such
    #or true positive rate           
    Sensitivity=100*(np.dot(complex_sizes.conj().T,Sn_co))/sum(complex_sizes)

    if (math.isnan(Sensitivity.tolist()[0])) or (math.isinf(Sensitivity.tolist()[0])):
        Sensitivity=0.001
    
    cluster_sum=sum(T)
    PPV=np.zeros((len(COMPLEXES_),len(new_clusters)))
    for i in range(len(COMPLEXES_)):
        PPV[i,:]=T[i,:]/cluster_sum
    for j in range(len(PPV)):
        for l in range(len(PPV[0])):
            if (math.isnan(PPV.tolist()[j][l])) or (math.isinf(PPV.tolist()[j][l])):
                PPV[j][l]=0

    PPV_cl=PPV.max(0) #cluster wise positive predictive value, #max item of each column

    #Positive predictive value (PPV, Precision) =Σ True positive/ Σ Test outcome positive
    Positive_Predictive_Value=100*(cluster_sum*np.matrix(PPV_cl).T)/sum(cluster_sum)
    if (math.isnan(float(Positive_Predictive_Value))) or (math.isinf(float(Positive_Predictive_Value))):
        Positive_Predictive_Value=0.001

    #Accuracy is the proximity of measurement results to the true value
    Accuracy=math.sqrt(Sensitivity*Positive_Predictive_Value) #geometrical mean

    #tropopoihsh kwdika g ypologismo seperation-pososto %100 opws k sta alla             

    sep_col=PPV
    sep_row=np.zeros((len(COMPLEXES_),len(new_clusters)))
    for j in range(len(new_clusters)):
        sep_row[:,j]=T[:,j]/np.sum(T, axis=1)
    for j in range(len(sep_row)):
        for l in range(len(sep_row[0])):
            if (math.isnan(sep_row.tolist()[j][l])) or (math.isinf(sep_row.tolist()[j][l])):
                sep_row[j][l]=0

    sep_matrix=sep_row*sep_col
    complex_wise_sep=np.sum(sep_matrix,axis=1)
    cluster_wise_sep=sum(sep_matrix)
    sep_co=sum(complex_wise_sep)/len(complex_wise_sep)
    sep_cl=sum(cluster_wise_sep)/len(cluster_wise_sep)
    Seperation=(math.sqrt(sep_co*sep_cl))*100


    return (Accuracy, Positive_Predictive_Value, Sensitivity, Seperation)
