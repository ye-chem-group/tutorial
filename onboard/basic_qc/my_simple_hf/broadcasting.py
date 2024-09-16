import sys
import numpy as np


if __name__ == '__main__':
    # A[i,j] = a[i] + b[j]  (len n)
    n = 3
    a = np.random.rand(n)
    b = np.random.rand(n)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = a[i] + b[j]
    A1 = a[:,None] + b
    print(A)
    print(A1)
