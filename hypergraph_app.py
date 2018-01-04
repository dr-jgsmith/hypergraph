#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 16:44:26 2017

@author: justinsmith
"""

from hypergraph import *

A = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

B = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

matrix = [[0, 1, 7, 0, 1, 6, 3, 2], 
          [5, 4, 7, 0, 0, 1, 4, 1],
          [1, 1, 1, 6, 6, 4, 3, 5],
          [1, 3, 7, 7, 3, 5, 5, 4],
          [1, 6, 6, 5, 4, 0, 0, 0],
          [1, 1, 1, 2, 4, 0, 0, 7],
          [0, 2, 2, 5, 0, 0, 0, 6],
          [0, 3, 3, 4, 4, 0, 0, 0],
          [0, 0, 7, 1, 0, 0, 1, 1],
          [2, 7, 2, 2, 0, 0, 1, 6],
          [1, 1, 1, 0, 4, 0, 0, 0],
          [3, 6, 6, 0, 1, 0, 0, 1]]


pre = [10, 12, 3, 10, 5, 19, 17, 0, 0, 9, 8, 1]
post = [7, 12, 3, 20, 8, 12, 17, 2, 1, 5, 6, 5]
w = [1.2, 0.5, 1, 1, 10, 5, 0.8, 1.2]



C = IncidentB(matrix, theta=3)
pd.DataFrame(C, index=A, columns=B)
simplex = get_simplex(C, A, B, 5)
print(simplex)
Conjugate = computeConjugate(C, A, B)
Conjugate[1]
faces = computeQFace(Conjugate[0], Conjugate[1], A, B)
QMatrix = computeQMatrix(faces, A)
cmps = computeQEqClass(faces, A)
print(faces)
print(QMatrix)
strct = computeQStruct(QMatrix)
print(strct)
chains = computeQChains(faces, strct, A)
print(chains)
ecc = Ecc(chains, strct)
print(ecc)
eccI = EccI(chains)
print(eccI)
comp = complexity(chains[1])
print(comp)
traffic = map_traffic(A, pre, post) #needs work and needs to be linked to q-transmission
print(traffic)
print(C)
psi = computePSI(C, w)
print(psi)
psin = computePSIN(psi)
print(psin)
norm = normPositive(matrix)
print(norm)
neg = normNeg(matrix)
print(neg)
f = computeQNear(QMatrix)
print(f)
#f2 = computeQNear2(Conjugate[1])