import numpy as np
import pickle
import math
#for 16 sites

def kit16():
    row = []
    col = []
    vals = []
    n = 16
    xlist=[[7],[10],[9],[4],[],[14],[13],[],[15],[],[],[12],[],[],[],[]]
    ylist=[[15],[2],[],[12],[11],[6],[],[8],[],[10],[],[],[],[14],[],[]]
    zlist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))

def kit16xy():
    row = []
    col = []
    vals = []
    n = 16
    xlist=[[7],[10],[9],[4],[],[14],[13],[],[15],[],[],[12],[],[],[],[]]
    ylist=[[15],[2],[],[12],[11],[6],[],[8],[],[10],[],[],[],[14],[],[]]
    zlist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[]]
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        for i in range(n):
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))

def kit16z():
    row = []
    col = []
    vals = []
    n = 16
    xlist=[[7],[10],[9],[4],[],[14],[13],[],[15],[],[],[12],[],[],[],[]]
    ylist=[[15],[2],[],[12],[11],[6],[],[8],[],[10],[],[],[],[14],[],[]]
    zlist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[]]
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a] += 1.0
                else:
                    vals[a] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
    return (vals,(row,col))

def kit162():
    row = []
    col = []
    vals = []
    n = 16
    xlist=[[7],[2],[],[4],[],[6],[],[],[9],[],[11],[],[13],[],[15],[]]
    ylist=[[1],[],[3],[],[5],[],[7],[],[15],[10],[],[12],[],[14],[],[]]
    zlist=[[8],[9],[10],[11],[12],[13],[14],[15],[],[],[],[],[],[],[],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))

def kit18():
    row = []
    col = []
    vals = []
    n = 18
    xlist=[[13],[6],[15],[8],[17],[10],[],[12],[],[14],[],[16],[],[],[],[],[],[]]
    ylist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[],[17],[]]
    zlist=[[5],[2],[],[4],[],[],[11],[8],[],[10],[],[],[17],[14],[],[16],[],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))
def kit24():
    row = []
    col = []
    vals = []
    n = 24
    xlist=[[11],[2],[],[4],[],[6],[],[8],[],[10],[],[],[13],[],[15],[],[17],[],[19],[],[21],[],[23],[]]
    ylist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[23],[14],[],[16],[],[18],[],[20],[],[22],[],[]]
    zlist=[[18],[13],[20],[15],[22],[17],[12],[19],[14],[21],[16],[23],[],[],[],[],[],[],[],[],[],[],[],[]]
    #to test
    #for i in range(n):
        #if xlist[i] != []:
            #print(i+1,xlist[i][0]+1)
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("kit24.pkl", "wb")
    #pickle.dump((vals,(row,col)),f)
    #f.close()
    return (vals,(row,col))

def kit242():
    row = []
    col = []
    vals = []
    n = 24
    xlist=[[11],[2],[],[4],[],[6],[],[8],[],[10],[],[],[13],[],[15],[],[17],[],[19],[],[21],[],[23],[]]
    ylist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[23],[14],[],[16],[],[18],[],[20],[],[22],[],[]]
    zlist=[[12],[13],[14],[15],[16],[17],[18],[19],[20],[21],[22],[23],[],[],[],[],[],[],[],[],[],[],[],[]]
    #to test
    #for i in range(n):
        #if xlist[i] != []:
            #print(i+1,xlist[i][0]+1)
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("kit242.pkl", "wb")
    #pickle.dump((vals,(row,col)),f)
    #f.close()
    return (vals,(row,col))

def kit243():
    row = []
    col = []
    vals = []
    n = 24
    xlist=[[11],[14],[13],[4],[],[18],[17],[8],[],[22],[21],[],[23],[],[],[16],[],[],[],[20],[],[],[],[]]
    ylist=[[23],[2],[],[16],[15],[6],[],[20],[19],[10],[],[12],[],[14],[],[],[],[18],[],[],[],[22],[],[]]
    zlist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[],[17],[],[19],[],[21],[],[23],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))

def kit244():
    row = []
    col = []
    vals = []
    n = 24
    xlist=[[7],[18],[9],[4],[],[22],[13],[],[15],[],[17],[12],[],[],[21],[],[23],[],[],[20],[],[],[],[]]
    ylist=[[23],[2],[],[12],[19],[6],[],[8],[],[10],[],[20],[],[14],[],[16],[],[18],[],[],[],[22],[],[]]
    zlist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[],[17],[],[19],[],[21],[],[23],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("indval.pkl", "wb")
    #pickle.dump((indices, vals),f)
    #f.close()
    return (vals,(row,col))

def kit245():
    row = []
    col = []
    vals = []
    n = 24
    xlist=[[17],[8],[19],[10],[21],[12],[23],[14],[],[16],[],[18],[],[20],[],[22],[],[],[],[],[],[],[],[]]
    ylist=[[1],[],[3],[],[5],[],[7],[],[9],[],[11],[],[13],[],[15],[],[17],[],[19],[],[21],[],[23],[]]
    zlist=[[7],[2],[],[4],[],[6],[],[],[15],[10],[],[12],[],[14],[],[],[23],[18],[],[20],[],[22],[],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**i+2**j)
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**i+2**j)
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    with open('kit245ham.pkl','wb') as f:
        pickle.dump((vals,(row,col)),f)
    #return (vals,(row,col))
#for 8
def kit8():
    row = []
    col = []
    vals = []
    n = 8
    xlist=[[3],[2],[],[],[5],[],[7],[]]
    ylist=[[1],[],[3],[],[7],[6],[],[]]
    zlist=[[4],[5],[6],[7],[],[],[],[]]
    old_off=0
    off=0
    for a in range(2**n):
        #if a % 10 ==0:
            #print(a)
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off=0
        for i in range(n):
            for j in zlist[i]:
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals[a+old_off] += 1.0
                else:
                    vals[a+old_off] -= 1.0
            #We get away with assigning the value to (a,b)
            #instead of (b,a) since the matrix is symmetrix
            for j in ylist[i]:
                b = a ^ (2**(i)+2**(j))
                if a_string[-(i+1)] == a_string[-(j+1)]:
                    vals.append(-1.0)
                else:
                    vals.append(1.0)
                row.append(a)
                col.append(b)
                off += 1
            for j in xlist[i]:
                b = a ^ (2**(i)+2**(j))
                vals.append(1.0)
                row.append(a)
                col.append(b)
                off +=1
    #f = open("kit8.pkl", "wb")
    #pickle.dump((vals,(row,col)),f)
    #f.close()
    return (vals,(row,col))


def mag001(n):
    row = []
    col = []
    vals = []
    for a in range(2**n):
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        for i in range(n):
            if a_string[-(i+1)]=='1':
                vals[a] += 1.0
            else:
                vals[a] -= 1.0
    return (vals, (row,col))

def mag111(n):
    row = []
    col = []
    vals = []
    old_off=0
    off=0
    for a in range(2**n):
        a_string = np.binary_repr(a, width = n)
        row.append(a)
        col.append(a)
        vals.append(0.0)
        old_off += off
        off = 0
        for i in range(n):
            #conjugating all elements (equivalent to transposing for hermitian matrixed)
            #since we're assigning values to (a,b) instead of (b,a)
            #(this doesn't actually matter, it's just to check)
            if a_string[-(i+1)]=='1':
                vals[a+old_off] += 1.0
                row.append(a)
                col.append(a ^ (2**i))
                vals.append(1.0-1.0j)
                off +=1
            else:
                vals[a+old_off] -= 1.0
                row.append(a)
                col.append(a ^ (2**i))
                vals.append(1.0+1.0j)
                off +=1
    #f = open('mag11124.pkl','wb')
    #pickle.dump((vals,(row,col)),f)
    #f.close()
    return (vals, (row,col))

#really -110
def mag110(n):
    row = []
    col = []
    vals = []
    for a in range(2**n):
        a_string = np.binary_repr(a, width = n)
        for i in range(n):
            #conjugating all elements (equivalent to transposing for hermitian matrixed)
            #since we're assigning values to (a,b) instead of (b,a)
            #(this doesn't actually matter, it's just to check)
            if a_string[-(i+1)]=='1':
                row.append(a)
                col.append(a ^ (2**i))
                vals.append(-1.0-1.0j)
            else:
                row.append(a)
                col.append(a ^ (2**i))
                vals.append(-1.0+1.0j)
    #f = open('mag11024.pkl','wb')
    #pickle.dump((vals,(row,col)),f)
    #f.close()
    return (vals, (row,col))

