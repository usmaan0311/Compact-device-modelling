# estimation of pi using monte carlo method
import numpy as np, sympy as sym
import matplotlib.pyplot as plt
import time


import random
inside = 0
i=1
n=int(input("Enter the total number of points: "))
it=time.time()
while (i<=n):
  x = random.random()
  y = random.random()
  if ((x**2)+(y**2))<=1:
    inside+=1
    plt.plot(x , y , 'go')
  else:
    plt.plot(x , y , 'ro')
  i+=1

ft=time.time()
pi=(4*inside)/n
print ("The value of pi is:")
print(pi)
print(f"The accuracy of estimated value of \u03C0 is :{(np.pi - pi)/np.pi}")
print(f"\nand it took {round(ft -it,4)} seconds to estimate it")
plt.show()
