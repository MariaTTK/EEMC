def buildMatrixW(NEIGHBORHOODS,NWEIGHTS,Weight):

    #Dhmiourgia araiou mhtrwou M (buildMatrixW))
    import numpy,time

    m=[[] for k in range(len(NWEIGHTS))]
    m1=[]

    for i in range(len(NWEIGHTS)): # sunolikos arithmos varwn gia kathe allhlepidrash
        m[i]=[i,len(NWEIGHTS[i][1][0])]
        m1.append([i]*m[i][1])

    import random 
    m1=[item for sublist in m1 for item in sublist]
    m1=[m1[0:int(len(m1)/2)],m1[int(len(m1)/2):int(len(m1))]]
    random.shuffle(m1[1])

    M=numpy.eye(len(NEIGHBORHOODS))

    for i in range(len(m1[0])):
       M[m1[0][i]][m1[1][i]]=Weight[i]
       M[m1[1][i]][m1[0][i]]=Weight[i]
    M=M/(numpy.ones(len(M))*sum(M))  
    M1=M

    return M1;

    '''#------------- PRINT A NUMPY MATRIX------------------------
    import numpy as np
    print('Print sparse matrix to file')
    t = time.time()

    # Write the array to disk
    with open("sparse_matrix_M1.txt", "wb") as f:
        # I'm writing a header here just for the sake of readability
        # Any line starting with "#" will be ignored by numpy.loadtxt
        f.write(b"#Array shape:{0}\n")

        # Iterating through a ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        for i in M1:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.  
            np.savetxt(f, i, fmt="%-7.4f")

            # Writing out a break to indicate different slices...
            f.write(b"   #New slice:\n")
            
    print(time.time() - t)
    #---------------------------------------------------------------'''
