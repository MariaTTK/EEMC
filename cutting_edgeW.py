def cutting_edgeW(INSIDEW,WDEGREES,clustersD,threshold):

    ratios=[]
    for i in range(len(clustersD)):
        s=[]
        for c in clustersD[i]:
            s.append(WDEGREES[c])
        outsideW=sum(s)
        try:
            ratios.append(INSIDEW[i]/(outsideW-INSIDEW[i]))
        except ZeroDivisionError:
            ratios.append(0)
    clustersCE=clustersD
    for r in range(len(ratios)):
        if (ratios[r]<threshold):
            clustersCE[r]=[]
    clustersCE= list(filter(None, clustersCE))

    return clustersCE
