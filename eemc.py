#global NEIGHBORHOODS WINDICES NWEIGHTS WDEGREES W COMPLEXES_
'''------------------------DATA-----------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------'''
    #Fortwsh arxeiou sumplokwn yeast
import time
import numpy as np
import math, sys

'''evaluation data'''
'''
file = open("BT_409.txt", "r") #opens file
complexes =(file.readlines()) #topothethsh sth lista complexes
file.close()
complexes_=[i.split() for i in complexes] #dhmiourgia listas COMPLEXES me sumploka, kathe sumploko mia lista me prwteines

    #diagrafh eggrafwn twn sumplokwn pou periexoun 1 prwteinh mono
for i in complexes_:
    if len(i)<2:
       del i[:]
COMPLEXES_=list(filter(None, complexes_))

    #afairesh duplicate grammwn tou dataset complexes
import itertools
COMPLEXES_.sort()
COMPLEXESS=list(COMPLEXES_ for COMPLEXES_,_ in itertools.groupby(COMPLEXES_))
COMPLEXES=list(filter(None,COMPLEXESS))
'''
inputfile=sys.argv[1]
#inputfile="gavin_yeast.eemc.input"
    #Fortwsh arxeiou YEAST with confidence scores
with open(inputfile, "r") as file: #opens file
        e =(file.readlines()) #topothethsh sth lista e
e_=[i.split() for i in e] #dhmiourgia listas e_ me sumploka kai confidence scores, kathe e_ mia lista me 2 prwteines kai scores

    #afairei an uparxei allhlepidrash metaksi ths idias prwteinhs
for i in range(len(e_)):
        if e_[i][0]==e_[i][2]:
                del e_[i]
                
    #afairesh duplicate grammwn tou dataset e_
import itertools
e_.sort()
E_=list(e_ for e_,_ in itertools.groupby(e_))
E=list(filter(None,E_)) #periexei se nest list ta "kathara" data
   

'''-------------------------NEIGHBOTHOODS, WEIGHTS------------------------------------
    ----------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------'''

Neighborhoods1, Neighborhoods2=[],[]
Interaction1, Interaction2=[],[]
Weights1, Weights2=[],[]
Igroups=[]
Weight=[]

for i in range(len(E)):
    Interaction1.append([E[i][0],E[i][2],E[i][3]]) #periexei tis allhlepidraseis twn 2 prwteinwn me vash th 1h sthlh
    Interaction2.append([E[i][2],E[i][0],E[i][3]]) #periexei tis allhlepidraseis twn 2 prwteinwn me vash th 2h sthlh
    Weight.append(float(E[i][3])) #periexei ta varh twn allhlepidrasewn
Interaction2.sort()    

    #Eyresh twn prwteinikwn allhlepidrasewn
import operator,itertools
from collections import defaultdict
from itertools import count

Ngroups1= itertools.groupby(Interaction1, key=operator.itemgetter(0)) 
Neighborhoods1=[[key,[y[1] for y in items]] for key,items in Ngroups1] #perixei listes:[prwteinh[prwteinh,..,prwteinh]]
Wgroups1= itertools.groupby(Interaction1, key=operator.itemgetter(0))
Weights1=[[key,list(y[2] for y in items)] for key,items in Wgroups1] #perixei listes:[prwteinh[varos,..,varos]]

Ngroups2= itertools.groupby(Interaction2, key=operator.itemgetter(0)) 
Neighborhoods2=[[key,[y[1] for y in items]] for key,items in Ngroups2] #perixei listes:[prwteinh[prwteinh,..,prwteinh]]
Wgroups2= itertools.groupby(Interaction2, key=operator.itemgetter(0))
Weights2=[[key,list(y[2] for y in items)] for key,items in Wgroups2] #perixei listes:[prwteinh[varos,..,varos]]

n=list(itertools.chain(Neighborhoods1, Neighborhoods2)) 
n.sort()
n1=itertools.groupby(n, key=operator.itemgetter(0))
n2=[[k,[y[1] for y in u]] for k, u in n1]

w=list(itertools.chain(Weights1, Weights2)) 
w.sort()
w1=itertools.groupby(w, key=operator.itemgetter(0))
w2=[[k,[y[1] for y in u]] for k, u in w1]

NEIGHBORHOODS=[[] for k in range(len(n2))]

   #Kanoume to NEIGHBORHOODS sth morphh ['a',['b',..,'d']]
for i in range(len(n2)):
   if len(n2[i][1])==2:
        NEIGHBORHOODS[i]=list(itertools.chain([n2[i][0],[n2[i][1][0]+n2[i][1][1]]]))
   else:
       NEIGHBORHOODS[i]=n2[i]

NWEIGHTS=[[] for k in range(len(w2))]

    #Kanoume to NWEIGHTS sth morphh ['a',['n1',..,'n3']]
for i in range(len(w2)):
   if len(w2[i][1])==2:
        NWEIGHTS[i]=list(itertools.chain([w2[i][0],[w2[i][1][0]+w2[i][1][1]]]))
   else:
       NWEIGHTS[i]=w2[i]

proteins=[]
for i in NEIGHBORHOODS:
	proteins.append(i[0])
	
WINDICES=[[] for i in range(len(NEIGHBORHOODS))]       
counts = defaultdict(lambda c=count(): next(c)) #The defaultdict() stores a new count() value each time a key has not yet been seen,
                                                #producing a unique counter per string.
for i in range(len(NEIGHBORHOODS)):
    seq=NEIGHBORHOODS[i][1]
    WINDICES[i]=list(itertools.chain([NEIGHBORHOODS[i][0],[list(a) for a in zip(*([counts[i] for i in sublist] for sublist in zip(*seq)))]]))

WDEGREES=[]
for i in range(len(WINDICES)):
    WDEGREES.append(sum([float(w) for w in NWEIGHTS[i][1][0]]))

'''evaluation data'''
'''
COMPLEXES_=[[[] for l in range(len(COMPLEXES[k]))] for k in range(len(COMPLEXES))]
s=[]

for i in COMPLEXES:
	s.append(len(i)) #periexei to megethos tou kathe sublist tou COMPLEXES

#linear time method:Get the indices of A where A ismember B. O(n+m) 
#for j in range(comp_size):
#   reverse_map = {x:i for i, x in enumerate(COMPLEXES[j])}
#    COMPLEXES_[j]=[(i) for i, x in enumerate(proteins) if x in reverse_map]
'''
ids=[[]for i in range(len(proteins))]
for i in range(len(proteins)):
	ids[i]=[proteins[i],counts[proteins[i]]] #to counts periexei ta keys(proteines) me ta antistoixa ids
'''evaluation data'''
'''
for idss in ids:
    for i in range(len(COMPLEXES)):
        for j in range(len(COMPLEXES[i])):
            if idss[0]==COMPLEXES[i][j]:
                COMPLEXES_[i][j]=idss[1]		
'''
#Dhmiourgia araiou mhtrwou M (buildMatrixW))
from buildMatrixW import buildMatrixW

M1=buildMatrixW(NEIGHBORHOODS,NWEIGHTS,Weight) #M1 5195x5195 sparse 
import scipy.sparse as sps
M1=sps.csr_matrix(M1)

last_generation=int(sys.argv[2])
chr_size=5
popul_size=int(sys.argv[3])
saved_r=-1
Pc=float(sys.argv[4])
print("Crossover parameter:", Pc)
Pm=float(sys.argv[5])
Fcon=int(sys.argv[6])
print("Adaptive Mutation Rate")
print("Mutation parameter:", Pm)
'''-------------------------XRWMOSWMATA------------------------------------
    ----------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------'''

initial_population=[[] for k in range(popul_size)]
A1,a1,a1_=[[] for k in range(popul_size)],[[] for k in range(popul_size)],[[] for k in range(popul_size)]
A2,a2,a2_=[[] for k in range(popul_size)],[[] for k in range(popul_size)],[[] for k in range(popul_size)]
A3,a3,a3_=[[] for k in range(popul_size)],[[] for k in range(popul_size)],[[] for k in range(popul_size)]
A4,a4,a4_=[[] for k in range(popul_size)],[[] for k in range(popul_size)],[[] for k in range(popul_size)]
A5,a5,a5_=[[] for k in range(popul_size)],[[] for k in range(popul_size)],[[] for k in range(popul_size)]
import numpy as np

for i in range(popul_size):
    initial_population[i]=np.random.choice([0, 1], size=(chr_size*10,)).tolist()
    

for i in range(popul_size): #SubLists
    a1[i]=initial_population[i][0:10]  #InflationRate(10bits)
    a2[i]=initial_population[i][11:21] #BestNeighborThreshold(10bits)
    a3[i]=initial_population[i][21:31] #HaircutParameter(10bits)
    a4[i]=initial_population[i][31:41] #DensityThreshold(10bits)
    a5[i]=initial_population[i][41:51] #CuttingEdgeThreshold(10bits)

    a1_[i]=[''.join((str(i) for i in a1[i]))] #%metatroph apo duadiko se mh arnhtiko dekadiko
    a2_[i]=[''.join((str(i) for i in a2[i]))]
    a3_[i]=[''.join((str(i) for i in a3[i]))]
    a4_[i]=[''.join((str(i) for i in a4[i]))]
    a5_[i]=[''.join((str(i) for i in a5[i]))]

    A1[i]=[int(i,base=2) for i in a1_[i]]
    A2[i]=[int(i,base=2) for i in a2_[i]]
    A3[i]=[int(i,base=2) for i in a3_[i]]
    A4[i]=[int(i,base=2) for i in a4_[i]]
    A5[i]=[int(i,base=2) for i in a5_[i]]
    
A=[A1, A2, A3, A4, A5]
A=[[[j/1024 for j in i] for i in l] for l in A]# 0<A<1

A[0]=[[j*2 for j in i] for i in A[0]] #InflationRate
A[4]=[[j*2 for j in i] for i in A[4]] #CuttingEdgeThreshold

'''
  A=[[A0,A0,..,A0] x10
     [A1,A1,..,A1]
     [A2,A2,..,A2]
     [A3,A3,..,A3]
     [A4,A4,..,A4]]

-------------------------MCL------------------------------------
    ----------------------------------------------------------------------------------
    -----------------------------------------------------------------------------------'''

fitness=np.zeros((popul_size,1))
fitness2=np.zeros((popul_size,1))
#to sizee_neighbour periexei ta mege8h ka8e geitonias
sizee_neighbour=[[] for k in range(len(WINDICES))]
for j in range(len(WINDICES)):
    sizee_neighbour[j]=[len(NEIGHBORHOODS[j][1][0])]#5195x2 sizee_neighbour, arithmos stoixeiwn kathe cell tou NEIGBORHOODS

'''Adaptive Mutation Rate'''
#Equation (9), dunamic control of the mutation parameters.
Pm_change=(Pm-(1/popul_size))/last_generation

best_elite=np.zeros((2,chr_size*10))
max_val_best=np.zeros((2,1))
dec_A=np.zeros((chr_size,1))
progress_elite=np.zeros((last_generation,1))
mean_progress_popul=np.zeros((last_generation,1))

print("Number of Generations:\n",last_generation)
print("Population size:\n", popul_size)

for generation in range(last_generation): #GENERATION
    #np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])
    from most_common import most_common
    groups=most_common(A[0])
    if saved_r!=groups:
        saved_r=groups
        saved_clusters=[]
    #find()
    a=np.matrix(A[0]).T
    common_r=np.nonzero(a==saved_r) #epistrefei p.x array([3, 6], dtype=int64)
                                    #A[0][3]==A[0][6]==saved_r

    for atomo in range(popul_size): #POPULATION SIZE
    #atomo=[0]
        atomo=[atomo]
    
        from MCL import MCL
        #ismember()
        if np.any(common_r==atomo[0]): #common_r is a list, epistrefei 1 pou brisketai h timh tou A sto B
            if (not saved_clusters)==True: #epistrefei 1 an to array einai keno
                for element in A[0][atomo[0]]:
                    r=float(element)+1.7
                clusters=MCL(NEIGHBORHOODS,M1,r)
                saved_clusters=clusters
            else:
                clusters=saved_clusters
        else:
            for element in A[0][atomo[0]]:
                    r=float(element)+1.7
            clusters=MCL(NEIGHBORHOODS,M1,r)
            
        #Filtering ( Good Neighbor, Hair Cut, Density, Cutting Edge)
        from filteringW import filteringW
        parameters=[]
        for p in A[1:(len(A))]:
            parameters.append(p[atomo[0]])
        print("Number of clusters before filtering:", len(clusters))
        new_clusters=filteringW(WINDICES, NWEIGHTS, WDEGREES,Weight,clusters,parameters)
        print("Perform filtering with Good Neighbor, Hair Cut, Density, Cutting Edge")
        print ("Number of clusters after filtering:", len(new_clusters))

        '''evaluation data'''
        '''
        from evaluation import evaluation
                Accuracy, Positive_Predictive_Value, Sensitivity, Seperation=evaluation(new_clusters,COMPLEXES_)
        '''
        if (not(new_clusters))==True:
            fitness[atomo[0]]=0.000000001
            continue

        '''-------------------------Genetic Algorithm------------------------------------
        ----------------------------------------------------------------------------------
        -----------------------------------------------------------------------------------'''
        #position,kathe lista einai o komvos kai kathe item h sustada pou ton periexei.
        #p.x position[5159][1]=398 --> new_clusters[398]=[.....,5159,....]

        #ALLHLEPIKALUPSH
        position=[[] for k in range(len(WINDICES))]

        for j in range(len(WINDICES)):
            for c in range(len(new_clusters)):
                    if np.any(np.array(new_clusters[c])==j):
                            position[j].append(c)     

        #EVALUATION
                            
        #Equation (5)
        #yv is calculated for each node v of cluster Cv
        posit_size=[]
        for j in range(len(WINDICES)):
            posit_size.append(len(position[j]))
        yv=np.zeros((len(WINDICES),1))
        yv_a=np.zeros((len(WINDICES), max(posit_size)))
        yv_b=np.zeros((len(WINDICES), max(posit_size)))
        yv_=[]
        for j in range(len(WINDICES)):
#periptwsh opou o komvos den anhkei se sustada    
            if not position[j]: #epistrefei True an to array einai keno
                yv[j]=sum([float(w) for w in NWEIGHTS[j][1][0]])

#periptwsh opou o komvos anhkei se sustada kai oi
#geitones tou den anhkoun sthn systada tou
            for l in range(posit_size[j]):
                    for k in range(len(WINDICES)):
                        if ((not position[j])==False): #to array DEN einai keno
                            if np.any(np.array(WINDICES[j][1][0])==k)==True:
                                if (np.any(np.array(new_clusters[position[j][l]])==k)==False): 
    
#yv_a:dianysma opou g ka8e komvo v periexei ta
#varh twn geitonwn ths kathgorias
                                   for gg in range(int(sizee_neighbour[j][0])):
                                       if  WINDICES[j][1][0][gg]==k:
                                          yv_a[j][l]=yv_a[j][l]+float(NWEIGHTS[j][1][0][gg])
                                

#periptwsh opou o komvos anhkei se sustada kai oi
#geitones tou anhkoun kai autoi sthn systada tou
                        if ((not position[j])==False): #to array DEN einai keno
                            if np.any(np.array(WINDICES[j][1][0])==k)==True:
                                if (np.any(np.array(new_clusters[position[j][l]])==k)==True):

#yv_b:dianysma opou g ka8e komvo v periexei ta
#1-varh twn geitonwn ths kathgorias
                                    for gg in range(int(sizee_neighbour[j][0])):
                                       if  WINDICES[j][1][0][gg]==k:
                                           yv_b[j][l]=yv_a[j][l]+(1-float(NWEIGHTS[j][1][0][gg]))

            yv_.append(np.sum([list(yv_a[j]),list(yv_b[j])],axis=0))

        yv_a=[]
        for j in range(len(WINDICES)):
            try:
                yv_a.append((sum(yv_[j]))/(len(np.nonzero(yv_[j]))))
            except ZeroDivisionError:
                yv_a.append(0)
        for i in range(len(yv_a)):
            if yv_a[i]!=0:
               yv[i]=np.array([yv_a[i]])
       
#Equation (6)
#bv is calculated for each node v of cluster Cv
        bv=np.zeros((len(WINDICES),1))
        bv_a=np.zeros((len(WINDICES), max(posit_size)))
        bv_b=np.zeros((len(WINDICES), max(posit_size)))
        bv_c=np.zeros((len(WINDICES), max(posit_size)))
        bv_=[]

        for j in range(len(WINDICES)):
            count_1,count_2,count_3=0,0,0
    #periptwsh opou o komvos den anhkei se sustada    
            if not position[j]: #epistrefei True an to array einai keno
                bv[j]=len(WINDICES[j][1][0])

#periptwsh opou o komvos anhkei se sustada kai oi
#geitones tou den anhkoun sthn systada tou
            for l in range(posit_size[j]):
                for k in range(len(WINDICES)):
                    if ((not position[j])==False): #to array DEN einai keno
                        if np.any(np.array(WINDICES[j][1][0])==k)==True:
                            if (np.any(np.array(new_clusters[position[j][l]])==k)==False):
                    #bv_a:periexei ton plh8os twn geitonwn tou
                    #komvou v ths kathgorias
                                count_1+=1
                                bv_a[j][l]=count_1

            #periptwsh opou o komvos anhkei se sustada kai oi
            #geitones tou anhkoun kai autoi sthn systada tou
                    if ((not position[j])==False): #to array DEN einai keno
                        if np.any(np.array(WINDICES[j][1][0])==k)==True:
                            if (np.any(np.array(new_clusters[position[j][l]])==k)==True):
                    #bv_b:periexei ton plh8os twn geitonwn tou
                    #komvou v ths kathgorias
                                count_2+=+1;
                                bv_b[j][l]=count_2

            #periptwsh opou o komvos anhkei se sustada kai oi
            #komvoi pou anhkoun sthn systada tou den einai
            #geitones tou
                    if ((not position[j])==False): #to array DEN einai keno
                        if np.any(np.array(WINDICES[j][1][0])==k)==False:
                            if (np.any(np.array(new_clusters[position[j][l]])==k)==True):
                    #bv_c:periexei ton plh8os twn geitonwn tou
                    #komvou v ths kathgorias
                                count_3+=+1;
                                bv_c[j][l]=count_3

    #ypologismos tou a8roismatos olwn twn komvwn (grammh) gia
    #oles tis systades stis opoies mporei na anhkoun (sthlh)
            bv_.append(np.sum([list(bv_a[j]),list(bv_b[j]),list(bv_c[j])],axis=0))
    
        bv_a=[]
        for j in range(len(WINDICES)):
            try:
                bv_a.append((sum(bv_[j]))/(len(np.nonzero(bv_[j]))))
            except ZeroDivisionError:
                bv_a.append(0)
        for i in range(len(bv_a)):
            if bv_a[i]!=0:
               bv[i]=np.array([bv_a[i]])

#Equation (7), minimazation function
        total_v=yv/bv
        fitness2[atomo[0]]=(((len(WINDICES)-1)/3)*(sum(total_v)))

#Equation (8), non-negative maximization function for the evolutionary algorithm
        fitness[atomo[0]]=(Fcon/(1+fitness2[atomo[0]])) #Fcon defined in line
        if (fitness[atomo[0]]<0):
            fitness[atomo[0]]=0
#end of range(popul_size)

    '''
    temp=np.setdiff1d(range(popul_size),range(len(A[0])) ) #Find the values in A that are not in B
    for item in temp:
        fitness[item]=fitness[item]
        warning('duplicate skipped')
    '''
    if generation!=0:
        min_index, min_value = min(enumerate(fitness), key=operator.itemgetter(1))
        if (max_value>min_value):
            initial_population[min_index]=elite
            fitness[min_index]=max_value
    max_index, max_value = max(enumerate(fitness), key=operator.itemgetter(1))
    elite=initial_population[max_index]

#elite: best_atomo trexousas genias
#eksetazetai to istoriko best elite atomo olwn twn geniwn kai
#apo8hkeuetai sto best_elite[1] me timh fitness sto
#max_val_best[1][0] kai decimal tou duadikou atomou sto dec_A
    best_elite[0]=elite
    max_val_best[0][0]=float(max_value[0])
    if max_val_best[0][0]>=max_val_best[1][0]:
        best_elite[1]=best_elite[0]
        max_val_best[1][0]=float(max_value[0])
        for i in range(chr_size):
            dec_A[i]=A[i][max_index]

    progress_elite[generation][0]=float(max_value[0])
    m=sum(fitness) / float(len(fitness))
    mean_progress_popul[generation][0]=m[0]


#convergence criterion-sugklish h sto last generation h sthn 5%
#sxetikh omoiothta apodoshs tou kalyterou atomou se sxesh me 
#thn mesh aposodh twn ypoloipwn atomwn

    sim_converge=abs(float(max_value[0])-m[0])/m[0]
    if (generation==last_generation) or (sim_converge<0.01):
        break

    oliki_apodosi=sum((np.array([(i/m) for i in fitness])**2))
    sxetiki_apodosi=(np.array([(i/m) for i in fitness])**2)/oliki_apodosi

    selected_indices=np.zeros((1,popul_size))
    for i in range(popul_size):
        temp=np.random.rand()
        for k in range(popul_size):
            if (temp<=sum(itertools.islice(sxetiki_apodosi, k))):
                selected_indices[0][i]=k
                break
#roulette wheel, afotou ginei elitismos, dld to
#xeirotero atomo antikatasta8ei me to kalutero
#
#tropopoihsh apo dekadiko se duadiko pali
#to B pleon periexei nai men plh8usmo tou A
#alla se duadikh morfh
#B:deksamenh zeugarwmatos

    "-------- Two Points Crossover Operator ---------"
    B=[]
    for i in selected_indices[0]:
            B.append(initial_population[int(i)])
	
#ypoxrewtikh diastaurwsh, an den epilegei kanena atomo gia
#diastaurwsh h diadikasia epanalamvanetai
    kk=0
    diast_atoma=[]
    diastaurwsh=[]
    while kk<=1:
        for i in range(popul_size):
            r=np.random.rand()
            if r<Pc:
                kk=kk+1
                diast_atoma.append(i)
                diastaurwsh.append(B[i]) #atoma pros diastaurwsh
            

#an epilegei perittos ari8mos, to teleutaio atomo aporriptetai
#afou h diastaurwsh ginetai se zeugaria
    if (kk%2)==1:
        del diastaurwsh[kk-1]
        kk=kk-1
    
#epilogh 2 shmeiwn diastaurwshs g ka8e atomo pros diastaurwsh
#ta shmeia den prepei na einai idia g ka8e sthlh k prepei na 
#exoun auksousa seira g ka8e grammh
    diastaurwshh=np.zeros((kk,50))
    shmeia_diastaurwshs=np.zeros((kk,2))
    def point(): #epistrefei, 1x2 array me tuxaia items kai sorted
        p=[round(elem[0]) for elem in (49*(np.random.rand(2,1))).tolist()]
        p.sort()
        return p
    i=0
    while (not(list(set(shmeia_diastaurwshs[:,0]).intersection(shmeia_diastaurwshs[:,1])))==False) and (i<len(shmeia_diastaurwshs)):
        p=point()
        shmeia_diastaurwshs[i]=p
        i=i+1   
    sizee=np.shape(diastaurwsh)
    for i in range(0,sizee[0],2):
        for j in range(int(shmeia_diastaurwshs[i,0])):
            diastaurwshh[i,j]=diastaurwsh[i][j]
            diastaurwshh[i+1,j]=diastaurwsh[i+1][j]
        for j in range(int(shmeia_diastaurwshs[i,0])+1,int(shmeia_diastaurwshs[i,1])):
            diastaurwshh[i,j]=diastaurwsh[i+1][j]
            diastaurwshh[i+1,j]=diastaurwsh[i][j]
        for j in range(int(shmeia_diastaurwshs[i,1])+1,sizee[1]):
            diastaurwshh[i,j]=diastaurwsh[i][j]
            diastaurwshh[i+1,j]=diastaurwsh[i+1][j]
    diastaurwsh=diastaurwshh

#sto B egine antikastash tou xeiroterou me
#to kalutero atomo k meta g ta atoma pou egine diastaurwsh
#egine antikastash autwn ston B
    for i in range(kk):
        B[diast_atoma[i]]=diastaurwsh[i]
    B=np.array(B)


    "-------- Mutation-----"
#ypologismos plh8ous omoiwn stoixeiwn metaksu tou elite kai twn
#allwn atomwn
    counter=np.zeros((popul_size,1))
    for i in range(popul_size):
        for j in range(50):
            if B[i][j]==elite[j]:
                counter[i]+=1
    counter=counter/50
    counter=((sum(counter))/(len(counter)))*100 #pososto % gia mesh omoiothta

    if counter>90:
        Pm+=Pm_change
        print("Mutation parameter:", Pm)
    else:
        Pm-=Pm_change
        print("Mutation parameter:", Pm)

    for i in range(popul_size):
        for j in range(50):
           r1=np.random.rand()
           if r1<Pm:
               if int(B[i][j])==1:
                   B[i][j]=0
               elif B[i][j]==0:
                    B[i][j]=1

#Metatroph tou B, apo numy array se list me integers                
    B=B.tolist()
    for i in range(len(B)):
            B[i] = [ int(x) for x in B[i] ]
    initial_population=B

#Metatroph tou A pali se dekadiko

    for i in range(popul_size): #SubLists
        a1[i]=initial_population[i][0:10]  #InflationRate(10bits)
        a2[i]=initial_population[i][11:21] #BestNeighborThreshold(10bits)
        a3[i]=initial_population[i][21:31] #HaircutParameter(10bits)
        a4[i]=initial_population[i][31:41] #DensityThreshold(10bits)
        a5[i]=initial_population[i][41:51] #CuttingEdgeThreshold(10bits)

        a1_[i]=[''.join((str(i) for i in a1[i]))] #%metatroph apo duadiko se mh arnhtiko dekadiko
        a2_[i]=[''.join((str(i) for i in a2[i]))]
        a3_[i]=[''.join((str(i) for i in a3[i]))]
        a4_[i]=[''.join((str(i) for i in a4[i]))]
        a5_[i]=[''.join((str(i) for i in a5[i]))]

        A1[i]=[int(i,base=2) for i in a1_[i]]
        A2[i]=[int(i,base=2) for i in a2_[i]]
        A3[i]=[int(i,base=2) for i in a3_[i]]
        A4[i]=[int(i,base=2) for i in a4_[i]]
        A5[i]=[int(i,base=2) for i in a5_[i]]
    
    A=[A1, A2, A3, A4, A5]
    A=[[[j/1024 for j in i] for i in l] for l in A]# 0<A<1

    A[0]=[[j*2 for j in i] for i in A[0]] #InflationRate
    A[4]=[[j*2 for j in i] for i in A[4]] #CuttingEdgeThreshold

#end for generation=1:last_generation

#paragwgh clusters apo to istoriko best elite olwn twn geniwn
r=dec_A[0][0]+1.7
clusters=MCL(NEIGHBORHOODS,M1,r)

parameters=[]
dec_A=dec_A.tolist()
for par in dec_A[1:(len(dec_A))]:
        parameters.append(par)
print("Number of clusters before filtering:", len(clusters))
new_clusters=filteringW(WINDICES, NWEIGHTS, WDEGREES,Weight,clusters,parameters)
print ("Number of clusters after filtering:", len(new_clusters))
'''evaluation data'''
'''
from evaluation import evaluation
Accuracy, Positive_Predictive_Value, Sensitivity, Seperation=evaluation(new_clusters,COMPLEXES_)
'''
#Metatroph twn id, pou uparxoun sto new_cluster, stis antistoixes proteines
clusters_NEW_proteins=[[[]for j in range(len(new_clusters[i]))]for i in range(len(new_clusters))]
for d in range(len(ids)):
    for i in range(len(new_clusters)):
        for j in range(len(new_clusters[i])):
            if ids[d][1]==new_clusters[i][j]:
                clusters_NEW_proteins[i][j]=ids[d][0]

outputfile=sys.argv[7]                
file=open(outputfile, "w")
file.writelines('Rank #Nodes Node labels\n')
for i in range(len(clusters_NEW_proteins)):
    n=len(clusters_NEW_proteins[i])
    file.writelines('%s\t'%i+'%s\t'%n+'\t'+'\t'+', '.join(clusters_NEW_proteins[i]).upper() + '\n')
file.close()
'''
noutputfile=sys.argv[2]
#outputfile="gavin_yeast.eemc.output"
cl=range(len(clusters_NEW_proteins))
with open(outputfile, "w") as file:
        for i in range(len(clusters_NEW_proteins)):
            for j in range(len(clusters_NEW_proteins[i])):
                file.writelines("cl_%s" %cl[i]+'\t'+clusters_NEW_proteins[i][j].upper()+'\n')
'''
