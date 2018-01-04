#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 08:20:42 2017

@author: justinsmith

"""

import numpy as np
import pandas as pd


def IncidentB(matrix, theta=1):
    """
    This function provides a basic method for describing a relation between two sets that
    have been computed as a MxN matrix of values that map to simplicies (rows) and vertices (columns).
    
    The theta value represents a threshold parameter for defining the partition of the matrix into
    0's and 1's.
    """
    B = np.zeros((len(matrix), len(matrix[0])))
    count = 0
    for i in matrix:
        cnt = 0
        for j in i:
            if j >= theta:
                B[count][cnt] = 1.0
            else:
                B[count][cnt] = 0.0
            cnt = cnt+1
        count = count+1
    return B


def get_simplex(B, row_headers, col_headers, theta):
    '''
    General function for looking displaying simplex level data.
    '''
    headers = ['simplex', 'num_ver', 'q_dim', 'slice'] + col_headers
    simplicies = []
    
    count = 0
    for i in B:
        tmp = []
        cnt = 0
        for j in i:
            if j == 1.:
                tmp.append(col_headers[cnt])
            else:
                pass  
            cnt = cnt + 1
        vertex = []
        for k in col_headers:
            if k in tmp:
                vertex.append(k)
            else:
                vertex.append('0')
        row = [row_headers[count], len(tmp), len(tmp)-1, theta] + vertex
        simplicies.append(row)
        count = count + 1
        
    dfSimplex = pd.DataFrame(simplicies, columns=headers)
    return dfSimplex


def computeConjugate(B, row_headers, col_headers):
    I = pd.DataFrame(B, index=row_headers, columns=col_headers)
    IT = I.T
    return I, IT

def computeQFace(I, IT, row_headers, col_headers):
    Z = I.dot(IT)
    return Z

def computeQMatrix(Z, row_headers):
    E = pd.DataFrame(np.ones(Z.shape), index=row_headers, columns=row_headers)
    S = Z.subtract(E)
    return S


def computeQEqClass(S, row_headers):
    '''
    This function
    
    '''
    cmp = []
    for i in S.values:
        tmp = []
        cnt = 0
        for j in i:
            if j > -1.0:
                tmp.append(row_headers[cnt])
            else:
                pass
            cnt = cnt + 1
        if len(tmp) == 0:
            pass
        else:
            cmp.append(tuple(tmp))
    x = sorted(set(cmp))
    new_x = [list(m) for m in x]
    test = []
    for m in new_x:
        d = []
        for n in m:
            for o in new_x:
                if n in o:
                    d = d + o
                else:
                    pass
        test.append(tuple(sorted(set(d))))
    return sorted(set(test))


def computeQStruct(S):
    matrix = np.array(S)
    QV = matrix.diagonal()
    return QV

def computeQChains(QFace, strct, row_headers):
    cmps = {}
    structure = {}
    q_strct = {}
    cnt = 0.0
    while max(strct) >= cnt:
        QFace = computeQMatrix(QFace, row_headers)
        cmp = computeQEqClass(QFace, row_headers)
        cmps[cnt] = cmp
        cnt = cnt + 1
    cnt = 0
    for j in row_headers:
        structure[j] = strct[cnt]
        cnt = cnt + 1
    for k in cmps.items():
        q_strct[k[0]] = len(k[1])
    return cmps, q_strct, structure


def EccI(q_chains):# eccI = 2(sum(q_dim/num_simps))/(q_dim*(q_dim+1))
    cmplx = {}
    for i in q_chains[2].items():
        n = int(i[1])
        occ = n+1
        x = list(range(occ))
        cnt = []
        for j in x:
            k = q_chains[0].get(float(j))
            for l in k:
                if i[0] in l:
                    d = len(l)
                else:
                    pass
            qi = j/d
            cnt.append(qi)
        
        if len(cnt) == 1:
            eccI = 0.0
        elif len(cnt) == 0:
            eccI = []
        else:
            eccI = (2*sum(cnt))/(max(x)*(max(x)+1))   
        cmplx[i[0]] = eccI
    return cmplx

    
def Ecc(q_chains, q_strct): #This is the mneasure of eccentricity presented by Casti
    cmplx = {}
    for i in q_chains[2].items():
        n = int(i[1])
        occ = n+1
        x = list(range(occ))
        data = {}
        for j in x:
            k = q_chains[0].get(float(j))
            for l in k:
                if i[0] in l:
                    d = len(l)
                else:
                    pass
            data[j] = d
        cmplx[i[0]] = data      
    ecc_casti = {}
    q_simplex = {}
    N = int(max(q_strct))+1
    f = list(range(N))
    count = 0
    for a in cmplx.items():
        tmp = 0
        for b in f:
            c = a[1].get(b)
            if c == None:
                pass
            else:
                if c > 1:
                    tmp = b
                else:
                    pass
        top_q = int(q_strct[count])
        bottom_q = tmp
        ecc = (top_q - bottom_q) / (bottom_q + 1)
        ecc_casti[a[0]] = ecc
        q_simplex[a[0]] = [top_q, bottom_q]
        count = count + 1
    return ecc_casti, q_simplex
           

def complexity(qvector):
    strct = []
    vect = []
    for i in qvector.items():
        x = i[0]+1
        y = x * i[1]
        strct.append(y)
        vect.append(i[1])
    z = sum(strct)
    compl = 2*(z/((max(vect)+1)*(max(vect)+2)))
    return compl


def map_traffic(row_headers, pre_vector_values, post_vector_values):
    data = {}
    count = 0
    for i in row_headers:
        x = pre_vector_values[count]
        y = post_vector_values[count]
        t_force = y-x
        data[i] = (x, y, t_force)
        count = count+1
    return data
            

def computePSI(q_simplex, weighted_vector):
    matrix = []
    for i in q_simplex:
        tmp = []
        cnt = 0
        for j in i:
            x = weighted_vector[cnt]
            d = x*j
            tmp.append(d)
            cnt = cnt + 1
        matrix.append(tmp)
    psi = np.array(matrix)
    return psi


def computePSIN(psi):
    psin = []
    for i in psi:
        x = max(i)
        row = []
        for j in i:
            if j == 0:
                p = 0
            else:
                if x == 0:
                    p = 0
                else:
                    p = j/x
            row.append(p)
        psin.append(row)
    psina = np.array(psin)
    return psina
            
    
#def computePCI():
    
    
#def computePCIN(pci):
    
    
def computeQNear(qmatrix):
    simplex = []
    sets = {}
    simplex_index = 0
    simplex = [i for i in qmatrix]
    for j in simplex:
        row = qmatrix.values[simplex_index]
        value_index = 0
        new_row = {}
        for l in row:
            if new_row.get(int(l)):
                new_row[int(l)].append(simplex[value_index])
            else:
                new_row[int(l)] = [simplex[value_index]]
            value_index =  value_index + 1
        X = [i[0] for i in new_row.items()]
        new_row2 = {}
        X = max(X)+1
        for i in list(range(X)):
            tmp = []
            for k in new_row.items():
                if i >= 0 and k[0] >= i:
                    [tmp.append(n) for n in k[1]]
                elif i == -1 and k[0] == i:
                    tmp = k[1]
                else:
                    pass
            new_row2[i] = tmp
        new_row2[-1] = new_row.get(-1)
        sets[j] = new_row2
        simplex_index = simplex_index + 1
    return sets
    

def computeQNear2(conjugate_face):
    simplex = []
    sets = {}
    X = dict(conjugate_face.sum())
    for i in dict(conjugate_face).items():
        print(i[0], dict(i[1]))
        

  
def normPositive(C):
    matrix = []
    for i in C:
        x = max(i)
        y = min(i)
        row = []
        for j in i:
            a = (j-y)/(x-y)
            row.append(a)
        matrix.append(row)      
    return matrix

def normNeg(C):
    matrix = []
    for i in C:
        x = max(i)
        y = min(i)
        row = []
        for j in i:
            a = (x-j)/(x-y)
            row.append(a)
        matrix.append(row)      
    return matrix
    
def get_thetas(matrix):
    nums = []
    for i in matrix:
        mn = min(i)
        nums.append(mn)
        mx = max(i)
        nums.append(mx)
    thetas = list(range(min(nums), max(nums)+1))
    return thetas
        
     
    
    
        
    
                
            
        
            
            
            
        
    
    
    
    
    
                    
            
        
        
        
        
        