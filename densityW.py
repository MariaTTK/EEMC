def densityW(WINDICES, Weight, clustersHC, threshold):
    import operator as op
    import functools
    density,INSIDEW=[],[]
    # nchoosek function,Binomial coefficient or all combinations n!/((n–k)! k!)
    def nchoosek(n, r):
        r = min(r, n-r)
        if r == 0: return 1
        numer = functools.reduce(op.mul, range(n, n-r, -1))
        denom = functools.reduce(op.mul, range(1, r+1))
        return numer//denom
    
    for i in range(len(clustersHC)):
        cluster=clustersHC[i]
        subset,subs,subs2,w=[],[],[],[]
        for c in cluster:
            subs.append(WINDICES[c][1][0])
        for lists in subs:
            for items in lists:
                subs2.append(items)
        #uniqify the list, order preserving method
        seen = {}
        for s in subs2:
            if s in seen: continue
            seen[s] = 1
            subset.append(s)
        for subset in subset:
            w.append(Weight[subset])
        INSIDEW.append(sum(w))

        nCk=nchoosek(len(cluster), 2)
    for insidew in range(len(INSIDEW)):
        density.append(INSIDEW[insidew]/nCk)
    clustersD=clustersHC
    for d in range(len(density)):
        if (density[d] <= threshold):
            INSIDEW[d]=[]
            clustersD[d]=[]
    clustersD= list(filter(None, clustersD))
    INSIDEW= list(filter(None, INSIDEW))
   
    
    return INSIDEW,clustersD
