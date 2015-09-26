def hair_cut2W(WINDICES, NWEIGHTS,clustersGN,param):
    #global WINDICES, NWEIGHTS
    import numpy as np
    from multiprocessing import Pool
    clustersHC=[[] for k in range(len(clustersGN))]
    
    for j in range(len(clustersGN)):
        cluster=clustersGN[j]
        connectivity=[]
        for c in range(len(cluster)):
            w=[]
            flags=[]
            for p in range(len(WINDICES[cluster[c]][1][0])):
                flags.append(np.any(np.array(cluster)==WINDICES[cluster[c]][1][0][p]))
            for f in range(len(flags)):
                if flags[f]==True:
                    w.append(NWEIGHTS[cluster[c]][1][0][f])
            connectivity.append(sum([float(w) for w in w]))       
            mean_conn=sum(connectivity) / float(len(connectivity))

        for i in range(len(connectivity)):
            if (connectivity[i] > (param*mean_conn)):
                clustersHC[j].append(cluster[i])

        #Infiltrations:Delete clusters with one item

        if (len(clustersHC[j])<2):
            clustersHC[j]=[]
    clustersHC= list(filter(None, clustersHC))

    return clustersHC
