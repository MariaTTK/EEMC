def good_neighborW(WINDICES, NWEIGHTS, WDEGREES,clusters,thres):
    
    #global WINDICES, NWEIGHTS, WDEGREES
    import numpy as np
    clustersGN=[[] for k in range(len(clusters))]    
    for j in range(len(clusters)):
        cluster = [list(row) for row in clusters[j]] #from tuple to a list
        
        neighbors=[]
        for c in cluster:
            neighbors.append(WINDICES[c[0]][1][0])
        neighbors=list(np.setdiff1d(neighbors[0], cluster[0])) #% returns the data in A that is not in B.


        flags=[]
        good_neighbors=[]
        adj=0
        for neighbor in neighbors:
            flags.append(any(p in WINDICES[neighbor][1][0] for p in cluster[0]))#ismember WINDICES,cluster
            for f in flags:
                if f==True:
                    adj=sum([float(w) for w in NWEIGHTS[neighbor][1][0]])

            if (adj/WDEGREES[neighbor] > thres):
                good_neighbors.append(neighbor)

        clustersGN[j]=cluster[0]+good_neighbors

    return clustersGN
