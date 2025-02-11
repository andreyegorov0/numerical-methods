import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 

def f(x,y):
    return np.array(np.e**x)

def Y1(x_n,y_n,h):
    k1 = h*f(x_n,y_n)
    k2 = h*f(x_n+2*h/3,y_n+2*k1/3)
    
    return y_n+(k1+3*k2)/4

def Y2(x_n,y_n,h):
    k1 = h*f(x_n,y_n)
    k2 = h*f(x_n + h/2,y_n + k1/2)
    k3 = h*f(x_n + 3*h/4,y_n + 3*k2/4)
    
    return y_n+(2*k1 + 3*k2 + 4*k3)/9 , Y1(x_n,y_n,h)

def Y3(x_n,y_n,h):
    k1 = h*f(x_n,y_n)
    k2 = h*f(x_n + h/2,y_n + k1/2)
    k3 = h*f(x_n + h/2,y_n + k2/2)
    k4 = h*f(x_n + h,y_n + k3)
    
    return y_n+(k1 + 2*k2 + 2*k3 + k4)/6 , 2*(k1 - k2 - k3 + k4)/3    

def cond(x_n,y_n,h,X,x_0,eps,method):
    if method==1:
        y1 = Y1(x_n,y_n,h)
        y2 = Y1(x_n+h/2, Y1(x_n,y_n,h/2) ,h/2)
        return np.max(np.abs(y2-y1)/(1-2**(-2)))<np.min(eps*h/(X-x_0))
    elif method==2:
        y2,y1 = Y2(x_n,y_n,h)
        return np.max(np.abs(y2-y1))<np.min(eps*h/(X-x_0))
    else:
        y1,E = Y3(x_n,y_n,h)
        return np.max(np.abs(E))<np.min(eps*h/(X-x_0))

def step(X,x_0,x_n,y_n,eps,method):
    h = X - x_n
    cnt = 1 
    while cond(x_n,y_n,h,X,x_0,eps,method)==False:
        h=h/2
        cnt+=1
    if method==1 or method==2:
        y_np1 = Y1(x_n,y_n,h)
    else:
        y_np1,_ = Y3(x_n,y_n,h)        
    return x_n+h, y_np1,cnt


def form(eps,x_0,X,y_0,method):
    x_n =  x_0
    y_n = y_0
    steps = 0
    while np.max(np.abs(x_n - X))>eps**2:
        if method==1:
            x_n, y_n,num = step(X,x_0,x_n,y_n,eps,method)
            steps+=6*num
        elif method==2:
            x_n, y_n,num = step(X,x_0,x_n,y_n,eps,method)
            steps+=5*num
        else:
            x_n, y_n,num = step(X,x_0,x_n,y_n,eps,method)
            steps+=4*num
    print('method №{0}:'.format(method))
    print('y(X) = {0}'.format(y_n))
    print('{0} calls of f(x,y)'.format(steps))


x_0 = np.array([0,0])
y_0 = np.array([1,1])
X = np.array([2,1])
eps = 10**(-4)
for i in range(1,4):
    form(eps,x_0,X,y_0,i)
