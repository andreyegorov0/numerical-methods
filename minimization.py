import numpy as np
import math
from numpy.linalg import norm
import matplotlib.pyplot as plt

eps = 10**(-8)
dim = 3
x_k = np.array([[5. , 5., 5.]],dtype='float64')

def f(x):
    return np.exp((x[0]-1)**2+10*(x[1]-2)**2+0.1*(x[2]-3)**2)


def Grad(x):
    tau = 0.1*math.sqrt(eps)
    gr = np.zeros(dim)
    for i in range(dim):
        x[i]+=tau
        f1 = f(x)
        x[i]-=2*tau
        f2 = f(x)
        x[i]+=tau
        gr[i] = (f1-f2)/(2*tau)
    return np.array(gr)


def check_1(part, x_k, x_km1):
    if part == 1:
        return norm(np.array(x_k)-np.array(x_km1), ord=2)<=math.sqrt(eps)
    else:
        return norm(np.array(x_k)-np.array(x_km1), ord=2)<=eps

def check_2(part, x_k, x_km1):
    if part == 1:
        return abs(f(x_k)-f(x_km1))<=math.sqrt(eps)
    else:
        return abs(f(x_k)-f(x_km1))<=eps

def check_3(part, x_k):
    if part == 1:
        return norm(Grad(x_k), ord=2)<=math.sqrt(eps)
    else:
        return norm(Grad(x_k), ord=2)<=eps

def H_1(x):
    return -Grad(x)

def Up(alpha,beta,x,h):
    mu=2
    prev=alpha
    alpha=mu*beta
    while f(x+alpha*h)<f(x+beta*h):
        prev=alpha
        alpha*=mu
    return prev

def Down(alpha,beta,x,h):
    lam=0.5
    alpha = beta*lam
    while f(x+alpha*h)>=f(x):
        alpha*=lam
    return alpha

def Alpha_1(h,x):
    beta = 0.1
    alpha = beta
    if f(x+alpha*h)<f(x):
        alpha = Up(alpha,beta,x,h)
    else:
        alpha = Down(alpha,beta,x,h)
    return alpha

def normalize(v):
    norma = norm(v,ord=2)
    if norma == 0: 
       return np.array(v)
    return np.array(v) / norma

def Alpha_2(h,x):
    n = 1
    a = np.array(x)
    b = np.array(x)+math.sqrt(eps)*normalize(h)
    a_n = np.array([a])
    b_n = np.array([b])
    x1_n = np.array([(a_n[-1]+b_n[-1])/2 - (eps/2**n)*normalize(h)])
    x2_n = np.array([(a_n[-1]+b_n[-1])/2 + (eps/2**n)*normalize(h)])
    while n<=np.log2((norm((b-a),ord=2)-eps)/eps):
        if f(x1_n[-1])<=f(x2_n[-1]):
            b_n = np.append(b_n,[x2_n[-1]],axis=0)
        else:
            a_n = np.append(a_n,[x1_n[-1]],axis=0)
        x1_n = np.append(x1_n, [(a_n[-1]+b_n[-1])/2 - (eps/2**n)*normalize(h)],axis=0)
        x2_n = np.append(x2_n, [(a_n[-1]+b_n[-1])/2 + (eps/2**n)*normalize(h)],axis=0)
        n+=1
    return (a_n[-1]+b_n[-1])/2

        

# def Second_der(x):
#     tau = 0.1*math.sqrt(eps)
#     g_xx = (f([x[0]+tau,x[1]]) +  f([x[0]-tau,x[1]]) - 2*f(x)) / (tau**2)
#     g_yy = (f([x[0],x[1]+tau]) + f([x[0],x[1]-tau]) - 2*f(x)) / (tau**2)
#     g_yx = (f([x[0]+tau,x[1]])-f(x) - f([x[0]+tau,x[1]-tau]) + f([x[0],x[1]-tau]))/ (tau**2)
#     gy_x = (f([x[0],x[1]+tau])-f(x) - f([x[0]-tau,x[1]+tau]) + f([x[0]-tau,x[1]]))/ (tau**2)
#     g_xy = 0.5 * (g_yx+gy_x) 
#     return np.array([[g_xx,g_xy],[g_xy,g_yy]])


def f_ij(i,j,x):
    tau = 0.1*math.sqrt(eps)
    x[j]+=tau
    f1=f(x)
    x[j]-=tau
    f2=f(x)
    x[j]+=tau
    x[i]-=tau
    f3=f(x)
    x[j]-=tau
    f4=f(x)
    x[i]+=tau
    return (f1-f2-f3+f4)/(tau**2)


def f_ii(i,x):
    tau = 0.1*math.sqrt(eps)
    x[i]+=tau
    f1 = f(x)
    x[i]-=2*tau
    f2 = f(x)
    x[i]+=tau
    f3 = f(x)
    num = (f1+f2-2*f3)/(tau**2)
    return num

def Gesse(x):
    matrix = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i,dim):
            if i==j:
                matrix[i][j]=f_ii(i,x)
            else:
                matrix[j][i] = 0.5*(f_ij(i,j,x)+f_ij(j,i,x))
                matrix[i][j] = matrix[j][i]
                
    return matrix

def H_2(x):
    return -(np.linalg.inv(Gesse(x))).dot(Grad(x))


cond=False
x_k = np.append(x_k,[x_k[-1]+10**(-12)],axis=0)

if (check_1(2,x_k[-1],x_k[-2]) == True) and (check_2(2,x_k[-1],x_k[-2]) == True) and (check_3(2,x_k[-1]) == True):
    print(x_k[-1])
else:
    t=0
    while cond==False:
        h = H_1(x_k[-1])
        part = 1
        t+=1
        alpha = Alpha_1(h,x_k[-1])
        x_k = np.append(x_k,[x_k[-1]+h*alpha],axis=0)
        if (check_1(part,x_k[-1],x_k[-2]) == True) and (check_2(part,x_k[-1],x_k[-2]) == True) and (check_3(part,x_k[-1]) == True):
            cond=True
            print('number of steps on the first half',t)
            print(x_k[-1])
    cond = False
    t=0
    while cond==False:
        part = 2
        t+=1
        h = H_2(x_k[-1])
        alpha = Alpha_2(h,x_k[-1])
        x_k = np.append(x_k,[alpha],axis=0)
        if (check_1(part,x_k[-1],x_k[-2]) == True) and (check_2(part,x_k[-1],x_k[-2]) == True) and (check_3(part,x_k[-1]) == True):
            print('final point = ',x_k[-1])
            print('number of steps on the second half',t)
            cond=True

# print(x_k)
# print(Grad(x_k[-1]))
# tau = 0.1*math.sqrt(eps)
# print((f([x_k[-1][0]+tau,x_k[-1][1]]) -  f([x_k[-1][0]-tau,x_k[-1][1]]) ))
