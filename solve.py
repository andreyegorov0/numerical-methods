import numpy as np

A  = np.loadtxt("A.txt", dtype=float)
B = np.loadtxt('B.txt',dtype=float)
C = np.linalg.inv(A) @ B
np.savetxt('X.txt', C, delimiter=' ')   