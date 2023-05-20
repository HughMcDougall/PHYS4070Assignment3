import numpy as np
import matplotlib.pyplot as plt

#------------------------
def kroneck(A,B):
    I, J = A.shape
    K, L = B.shape

    out = np.zeros([I*K, J*L])

    for i in range(I):
        for j in range(J):

            for k in range(K):
                for l in range(L):

                    out[i*K + k , j*L + l] = A[i,j] * B[k,l]

    return(out)

def kp(X):
    N = len(X)
    out = X[N-1]
    for i in range(1,N):
        out = kroneck(X[N-1-i], out)
    return(out)

#------------------------
N = 3

I = np.eye(2)
A = np.zeros([2,2])
B = np.zeros([2,2])
A[0,0], A[1,1] = 1/2, -1/2
B[0,1], B[1,0] = 1/2,  1/2

def sigma_z(m):
    X = [np.eye(2**(m)), A, np.eye(2**(N-1-m))]
    return( kp(X) )

def sigma_x(m):
    if m==N: return(sigma_x(0))
    X = [np.eye(2**(m)), B, np.eye(2**(N-1-m))]
    return( kp(X) )

def H(g=0):
    
    out = np.zeros([2**N,2**N])

    for m in range(N):
        out-= sigma_z(m)
        out+= - g * np.dot(sigma_x(m),sigma_x(m+1))
    return(out)


#--------------------
fig,ax = plt.subplots(3,1)
ax[0].imshow(H(g=0))
ax[1].imshow(H(g=1))
ax[2].imshow(H(g=2))
plt.show()
