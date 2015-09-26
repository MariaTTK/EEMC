def MCL(NEIGHBORHOODS,M1,r):
    import numpy as np
    import scipy.sparse as sps
    import time,struct
    epoch=0
    change=True
    buffer=np.empty([len(NEIGHBORHOODS),len(NEIGHBORHOODS)])
    Mdensity=[]
    

    while (change and epoch<150): #argei
        #t=time.time() 
        epoch=epoch+1
        if (epoch>25):
            if (np.all(M1 == buffer)):
                import warnings
                warnings.warn("oscillation detected")
                break
            else:
                buffer=M1
        #t1=time.time()        
        M2=M1*M1 # Ektonosi, pol/mos tou M me ton eauto tou
        #print("M2=M1*M1",time.time() - t1)
        Mdensity.append(len(np.nonzero(M2)[1])/np.size(M2)) #Number of nonzero matrix elements(M2)/Number of array elements
        if Mdensity[epoch-1]>0.2:
            M2=M2.todense()
            M2[M2<0.00025]=0
            Mdensity[epoch-1]=len(np.nonzero(M2)[1])/np.size(M2) #nonzero,find the indices where a!=0
            if (Mdensity[epoch-1]<=0.2):
                M2=sps.csr_matrix(M2)

            S=sum(M2) #sum kathe stilis tou M2

            M2=M2/(np.ones((np.size(S),1))*S)        
        #t2=time.time()
        from inflation import inflation
        M1=inflation(M2,r)
        #print("inflation(M2,r)",time.time() - t2)
        c=np.not_equal(M1,M2)
        change=np.any(c)
    #print("While",time.time() - t)
    #t3=time.time()
    for i in range(1,len(NEIGHBORHOODS)+1):
            if (list(np.shape(np.nonzero(M1[:,i-1])))[1])>0:
                    indices=np.nonzero(M1[:,i-1])
                    z=np.zeros(((len(indices[0]))-1 , 1 ))
                    uv=[1]
                    c=-1
                    for j in z:
                            uv.append(z[0][0])
                    for l in indices[0]:
                            c=c+1
                            M1[l,i-1]=uv[c] 
    #print("M1 me asous",time.time() - t3)
    flags=np.ones((1,len(NEIGHBORHOODS)), dtype=bool)  #array of logical ones.
    clusters=[[] for i in range(len(NEIGHBORHOODS))]
    t4=time.time()
    for i in range(1,(len(NEIGHBORHOODS)+1)):
        if M1[i-1,i-1]: #if M1[i-1,i-1]!=0
            uint16=struct.unpack('H', struct.pack('h',i-1))
            clusters[i-1].append(uint16)
            flags[0][i-1]=0
    #print("1o clusters",time.time() - t1)

    pick=np.array(range(1,len(NEIGHBORHOODS)+1))
    t5=time.time()
    for i in pick[flags[0]]:
        for j in range(1,len(clusters)+1):
            if np.any(M1[clusters[j-1],i-1]):
                uint16=struct.unpack('H', struct.pack('h',i-1))
                clusters[j-1].append(uint16)
                flags[0][i-1]=0
                break
    #print("2o clusters",time.time() - t5)
    clusters=list(filter(None, clusters))

    return clusters
