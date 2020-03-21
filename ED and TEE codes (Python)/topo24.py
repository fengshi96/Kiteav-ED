import numpy as np
import numpy.linalg
import kit_ham
import scipy.sparse
import scipy.sparse.linalg as lg
import math
import time
import pickle
import sys

i = sys.argv[1]
num_points = 100
minrange = 0
maxrange = 1.57
step = (maxrange-minrange)/ num_points
n = 24
dim = 2**n
lengthA = 6
lengthB = lengthA
lengthC = lengthA
#could just do A+B
lengthAB = 12
lengthBC = lengthAB
lengthAC = lengthAB
lengthABC = 6

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
    return split(full_list,[6,7,8,12,13,14])

def splitB(full_list):
    return split(full_list,[9,10,11,21,22,23])

def splitC(full_list):
    return split(full_list,[0,1,2,18,19,20])

def splitAB(full_list):
    return split(full_list,[6,7,8,9,10,11,12,13,14,21,22,23])

def splitBC(full_list):
    return split(full_list,[0,1,2,9,10,11,18,19,20,21,22,23])

def splitAC(full_list):
    return split(full_list,[0,1,2,6,7,8,12,13,14,18,19,20])

def splitABC(full_list):
    return split(full_list,[3,4,5,15,16,17])

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

with open('/home/vengal.8/kit24/kit24mag111' + i +'.pkl','rb') as t:
    weights = pickle.load(t)
SA = ee(reduced_eigs(rweights(weights, splitA, lengthA),lengthA))
SB = ee(reduced_eigs(rweights(weights, splitB, lengthB),lengthB))
SC = ee(reduced_eigs(rweights(weights, splitC, lengthC),lengthC))
SAB = ee(reduced_eigs(rweights(weights, splitAB, lengthAB),lengthAB))
SBC = ee(reduced_eigs(rweights(weights, splitBC, lengthBC),lengthBC))
SAC = ee(reduced_eigs(rweights(weights, splitAC, lengthAC),lengthAC))
SABC = ee(reduced_eigs(rweights(weights, splitABC, lengthABC),lengthABC))
with open('kit24mag111.txt','a') as t:
    print(i,SA+SB+SC-SAB-SBC-SAC+SABC,file=t)




