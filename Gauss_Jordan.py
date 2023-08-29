import numpy as np 
import cmath
import math
from numpy import linalg as LA

tol = 0.00000000000001

def get_matrix(): #Doc ma tran A tu file
	f = open('Gauss_Jordan.txt')
	tmp_A = [[], [], [], [], [], [], [], [], [], []]
	tmp_b = [[], [], [], [], [], [], [], [], [], []]
	i = 0
	j = 0
	for line in f:
		number_array = line.split()
		if len(number_array) == 2:
			SoHang = float(number_array[0])
			SoCot = float(number_array[1])
		elif len(number_array) == 1:
			tmp_b[j].append(float(number_array[0]))
			j += 1
		else:
			for number in number_array:
				tmp_A[i].append(float(number))
			i += 1
	A = np.array(tmp_A)
	b = np.array(tmp_b)
	return A, b

def gj(A, b, tol):
    m, n = A.shape
    i = 0
    j = 0
    A = np.concatenate((A, b), axis=1)
    while i < m and j < n:
        # Find value and index of largest element in the remainder of column j.
        k = np.argmax(np.abs(A[i:, j])) + i
        p = np.abs(A[k, j])
        print("****")
        print(A[:, n])
        print("****") 
        if p < tol:
            # The column is negligible, zero it out.
            A[i:, j] = 0
            j += 1
        else:
            # Swap i-th and k-th rows.      
            A[[i, k], j:n+1] = A[[k, i], j:n+1]
            # Divide the pivot row by the pivot element.
            A[i, j:n+1] = A[i, j:n+1] / A[i, j]
            # Subtract multiples of the pivot row from all the other rows.
            for k in np.concatenate((np.arange(i), np.arange(i+1, m))):
                A[k, j:n+1] = A[k, j:n+1] - A[k, j] * A[i, j:n+1]
            i += 1
            j += 1
            print("****")
            print(A[:, n])
            print("****")      
    return A

A, b = get_matrix()
E = np.eye(len(A))
A = E - A
print(A)
print(b)
print(LA.inv(A)@b)

print()
M = gj(A, b, tol)
print(M)		