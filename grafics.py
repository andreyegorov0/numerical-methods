import numpy as np
import matplotlib.pyplot as plt
import math
 
filefile = 'data.txt'

with open(filefile, 'r') as file:
    file.seek(0)
    arrXi = [float(x) for x in file.readline().split()]
    arrYi = [float(x) for x in file.readline().split()]
    arrNx = [float(x) for x in file.readline().split()]
    arrNy = [float(x) for x in file.readline().split()]
 

# print("first arrXCCCCCCCi:", arrNx)
# print("second arrYi:", arrNy)
plt.title("Interpolation polynom")
#plt.scatter(arrNx, arrNy)
plt.plot(arrNx, arrNy, label = 'interp',color='orange' ) #первое массив иксов, второе массив игриков
plt.scatter(arrXi, arrYi)
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
p = np.arange(-6.5,6.51,0.01) #6.5 заменил на 6.51, в функции range левая граница включается а правая нет + разбиение помельче сделал 0.1 -> 0.01
y = np.array([math.sin(a) for a in p])
plt.plot(p,y,label = 'func', color='red')
plt.legend(loc='best')
plt.xlabel("x")
plt.ylabel("y")
plt.show()

  

 # нельзя через нью нужно маллок нельзя через стек нужно основеуб память