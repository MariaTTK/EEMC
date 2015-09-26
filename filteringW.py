def filteringW(WINDICES, NWEIGHTS, WDEGREES,Weight,clusters,parameters):
    
    #global WINDICES, NWEIGHTS, WDEGREES

    from good_neighborW import good_neighborW
    clustersGN=good_neighborW(WINDICES, NWEIGHTS, WDEGREES,clusters,parameters[0][0])

    from hair_cut2W import hair_cut2W
    clustersHC=hair_cut2W(WINDICES, NWEIGHTS,clustersGN,parameters[1][0])

    from densityW import densityW
    INSIDEW,clustersD=densityW(WINDICES,Weight,clustersHC,parameters[2][0])

    from cutting_edgeW import cutting_edgeW
    clustersCE=cutting_edgeW(INSIDEW,WDEGREES,clustersD,parameters[3][0])
    
    new_clusters=clustersCE
    return new_clusters
   
