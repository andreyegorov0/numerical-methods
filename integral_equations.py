import numpy as np
import math
import pylab
from matplotlib import mlab
import matplotlib.pyplot as plt


eps = 10**(-3)
a, b = 0. , 1.
# a, b = 0.01, 0.5 #2
# a, b = 1e-6, 2*math.pi #3/5
n = 10
PLOT_N = 100



def exact_solution(x):
    return x + np.exp(-x) #1
    # return 1.0 #2
    # return (25.0 + 27.0*np.cos (2*x)) / (160*math.pi) #3
    # return np.cos (2*x) #5
    # return x #8


def K(x, t):
    return x*np.exp(t)/2 #1
    # return np.sin (x*t) #2
    # return -1 / (4 * math.pi * (np.sin ((x + t) / 2)**2 + 0.25 * np.cos ((x + t)/2)**2)) #3
    # return -np.sin (x)*np.cos (t) #5
    # return x*t/2


def f(x):
    return np.exp(-x) #1
    # return 1.0 + (np.cos (x/2) - 1)/x #2
    # return (5.0 + 3.0*np.cos (2*x)) / (16*math.pi) #3
    # return np.cos (2*x) #5
    # return 5*x/6


def mesh(a, b, n):
    h_mesh = (b-a) / n
    return np.array([i * h_mesh for i in range(n+1)], dtype='float64')

def integral(x, y, g, h_1):
    return (g(x) + g(y)) * h_1 / 2

def integral_2(x, y, g, h_1):
    return (g(x) + 2*g((x+y)/2) + g(y)) * h_1 / 4

def Step(x, y, a, b, g, h_1):
    x_1, y_1 = x, y
    h_2 = h_1
    I_h = integral(x_1, y_1, g, h_1)
    I_h2 = integral_2(x_1, y_1, g, h_1)
    if np.abs(I_h-I_h2) <= eps * h_1 / (2*(b-a)):
        return h_2, I_h2
    else:
        cnt = 0
        while (np.abs(I_h - I_h2) > (eps * h_2 / (2*(b-a))) and cnt < 30):
            h_2 /= 2
            y_1 = x_1 + h_2
            cnt += 1
            I_h = integral(x_1, y_1, g, h_2)
            I_h2 = integral_2(x_1, y_1, g, h_2)
        return h_2, I_h2


def l2_norm(a, b, g): 
    N_1 = 20
    h_1 = (b-a) / N_1
    S_1 = 0
    x_2 = a 
    y_2 = x_2 + h_1
    while np.abs(x_2-b) > eps**2:
        step, s = Step(x_2, y_2, a, b, g, h_1)
        S_1 += s
        x_2 += step
        if (x_2 + h_1) < b:
            y_2 = x_2 + h_1
        else:
            y_2 = b 
            h_1 = y_2 - x_2
    return np.sqrt(S_1)


def A(a, b, n): #матрица коэффициентов
    A_1 = np.zeros(n+1)
    h_a = (b-a) / (n)
    for i in range(n+1):
        if i == 0 or i == n:
            A_1[i] = h_a / 2
        else:
            A_1[i]=h_a
    return np.array(A_1, dtype=np.float64)

def K_i(x, a, b, n): #разбиение ядра на узлах
    vector_ki = []
    m_ki = mesh(a, b, n)
    for i in range(n+1):
        y = m_ki[i]
        k_ki = K(x,y)
        vector_ki.append(k_ki)
    return  np.array(vector_ki, dtype=np.float64)
    
def f_i(a, b, n):#Разбиение функии f на узлах 
    vector_fi = []
    m_fi= mesh(a, b, n)
    for i in range(n+1):
        y = m_fi[i]
        k_fi = f(y)
        vector_fi.append(k_fi)
    return  np.array(vector_fi, dtype=np.float64)


def K_matrix(a, b, n):
    k_matrix = np.zeros((n+1,n+1))
    m_km = mesh(a,b,n)
    for i in range(n+1):
        for j in range(n+1):
            k_matrix[i][j] = K(m_km[j],m_km[i])
    return k_matrix

def u_n(a, b, n): #возвращает решение слау
    fn = f_i(a,b,n)
    n_un = n
    matrix_k_un = K_matrix(a, b, n).T
    h_un = (b-a) / n_un
    for i in range(n_un+1):
        for j in range(n_un+1):
            if j == 0 or j == n_un:
                matrix_k_un[i][j] *= h_un/2
            else:
                matrix_k_un[i][j] *= h_un
    U_n = np.zeros(n_un+1)
    U_n = (np.linalg.inv(np.eye(n_un+1) - matrix_k_un)).dot(fn)
    return np.array(U_n,dtype='float64')

def u(x,a,b,u_1): #считает u в произвольной точке
    N = u_1.shape [0]
    c = A(a, b, N)
    k_u = K_i(x, a, b, N)
    S_u = 0
    for i in range(N):
        s_u = c[i]*k_u[i]*u_1[i]
        S_u += s_u
    result_u = f(x) + S_u
    return result_u





u_1 = u_n(a, b, n)
n *= 2
u_2= u_n(a, b, n)
print('L_2 norm = ', l2_norm(a, b, lambda x: (u(x, a, b, u_1) - u(x, a, b, u_2))**2))
print('Узлов: {0} '.format(n+1))



while l2_norm(a, b, lambda x: (u(x, a, b, u_1) - u(x, a, b, u_2))**2)>eps:
    u_1 = u_2
    n *= 2
    u_2 = u_n(a, b, n)
    print('L_2 norm = ', l2_norm(a, b, lambda x: (u(x, a, b, u_1) - u(x, a, b, u_2))**2))
    print('Узлов: {0} '.format(n+1))

print ("Отклонение по норме в L2 от точного решения:")
print (  l2_norm (a, b, lambda x: (u(x, a, b, u_2) -  exact_solution (x))**2))

xList = mesh(a, b, PLOT_N)
uList = np.array ([u(x, a, b, u_1) for x in xList])
sList = np.array ([exact_solution (x) for x in xList])
devList = np.abs (uList - sList)
maxDev = devList.max ()
    
print ("Отклонение по норме в C [{0}, {1}] от точного решения:".format (a, b))
print (maxDev)
pylab.plot (xList, uList, color = 'red')
pylab.plot (xList, sList, color = 'black', linestyle = 'dashed')
plt.show()
