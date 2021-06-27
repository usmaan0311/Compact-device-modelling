from pi_MonteCarlo_parallel import*


n=input("enter number of points:\n")
n=int(n)

it=time.time()
p=monte_carlo_pi_part(n)
ft=time.time()
pii=p*4/n

print(f"Esitmated value of Pi: {pii} ")  
print(f"Accuracy of estimated value of \u03C0 is : {(pii - math.pi)/math.pi} \n and it took {round(ft -it,4)} seconds to estimate it")
