import numpy as np
import numpy.linalg
import kit_ham
import scipy.sparse
import scipy.sparse.linalg as lg
import math
import time
import pickle
import sys

i = int(sys.argv[1])
num_points = 500
minrange = 0
maxrange = 1.57
step = (maxrange-minrange)/ num_points
n = 18
dim = 2**n
lengthA = 5
lengthB = 4
lengthC = lengthB
#could just do A+B
lengthAB = 9
lengthBC = 8
lengthAC = lengthAB
lengthABC = 5

def frange(start,stop,step):
    i = start
    while i< stop:
        yield i
        i += step

def split(full_list, nums):
    list1 = full_list[nums[0]]
    for i in nums[1:]:
        list1 += full_list[i]
    list2 = full_list[:nums[0]]
    length = len(nums)
    for i in range(length-1):
        list2 += full_list[(nums[i]+1):nums[i+1]]
    list2 += full_list[nums[-1]+1:]
    return [list1, list2]

def splitA(full_list):
    return split(full_list,[0,1,2,3,4])

def splitB(full_list):
    return split(full_list,[6,7,8,9])

def splitC(full_list):
    return split(full_list,[12,13,14,15])

def splitAB(full_list):
    return split(full_list,[0,1,2,3,4,6,7,8,9])

def splitBC(full_list):
    return split(full_list,[6,7,8,9,12,13,14,15])

def splitAC(full_list):
    return split(full_list,[0,1,2,3,4,12,13,14,15])

def splitABC(full_list):
    return split(full_list,[5,10,11,16,17])

#alength is number of sites in reduced region
def rweights(weights,split,alength):
    bdim = 2**(n-alength)
    rweights = [[] for _ in range(bdim)]
    for w in weights:
        splita,splitb = split(w[1])
        rweights[int(splitb,2)].append([w[0],splita])
    return rweights

def reduced_eigs(rweights,alength):
    dim = 2**alength
    rho = np.zeros((dim,dim),dtype=complex)
    for rweight in rweights:
        for w in rweight:
            for x in rweight:
                rho[int(w[1],2),int(x[1],2)]+=w[0]*np.conj(x[0])
    return numpy.linalg.eigvalsh(rho)

def ee(eigs):
    lg = np.vectorize(lambda x: 0 if x<1e-6 else math.log2(x),otypes =[np.float])
    return -np.dot(eigs,lg(eigs))

K = scipy.sparse.csr_matrix(kit_ham.kit18(),shape=(dim,dim))
h = scipy.sparse.csr_matrix(kit_ham.mag111(18),shape=(dim,dim))
for r in range(10*i,10*(i+1)):
    q = r * step
    v = lg.eigsh(np.cos(q)*K-np.sin(q)/math.sqrt(3)*h,1,which='SA')[1]
    v = v[:,0]/numpy.linalg.norm(v[:,0])
    weights = [[v[i],np.binary_repr(i,n)] for i in range(dim) if numpy.linalg.norm(v[i]) > 1e-4]
    SA = ee(reduced_eigs(rweights(weights, splitA, lengthA),lengthA))
    SB = ee(reduced_eigs(rweights(weights, splitB, lengthB),lengthB))
    SC = ee(reduced_eigs(rweights(weights, splitC, lengthC),lengthC))
    SAB = ee(reduced_eigs(rweights(weights, splitAB, lengthAB),lengthAB))
    SBC = ee(reduced_eigs(rweights(weights, splitBC, lengthBC),lengthBC))
    SAC = ee(reduced_eigs(rweights(weights, splitAC, lengthAC),lengthAC))
    SABC = ee(reduced_eigs(rweights(weights, splitABC, lengthABC),lengthABC))
    with open('/home/vengal.8/kit18mag111p1.txt','a') as t:
        print(r,SA+SB+SC-SAB-SBC-SAC+SABC,file=t)


"""pd.DataFrame(data).to_csv('kit8topo111.csv', header=None)
print(time.time()-time1)"""




