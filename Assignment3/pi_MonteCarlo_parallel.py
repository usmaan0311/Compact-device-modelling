import math
import random
import multiprocessing
import time
from multiprocessing import Pool



def monte_carlo_pi_part(n):
    
    count = 0
    for i in range(n):
        x=random.random()
        y=random.random()
        
        
        if x*x + y*y <= 1:
            count=count+1
        
    
    return count


if __name__=='__main__':
    
    np = multiprocessing.cpu_count()
    

    
    n = input("Enter the number of points:\n")
    n=int(n)
    
   
    part_count=[n//np for i in range(np)]

    it=time.time()
   
    pool = Pool(processes=np)   

    ft=time.time()    
  
    count=pool.map(monte_carlo_pi_part, part_count)
    pii=sum(count)/(n*1.0)*4

    print( "Esitmated value of Pi:: ", sum(count)/(n*1.0)*4)  
    print(f"Accuracy of estimated value of \u03C0 is : {(pii - math.pi)/math.pi} \n and it took {round(ft -it,4)} seconds to estimate it") 
