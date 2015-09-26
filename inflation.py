def inflation(M,r):
    import numpy as np
    temp=np.power(M,r) #shape(temp)=(5195, 5195)

    #allagh gia dior8wsh out of memory
    S=sum(temp) #shape(S)=(1, 5195)


    G=temp/(np.ones((np.size(S),1))*S)

    M1=G

    return M1;
