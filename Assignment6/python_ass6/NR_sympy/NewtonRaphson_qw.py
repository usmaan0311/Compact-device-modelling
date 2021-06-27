import sympy as sym

#guess function: f=3*x**2 -4*(1-x)*exp(x) -sin(x)*cos(x)
#s=input("enter function of your choice in variable x\n")









# Newton Raphson method
#s=" x**2 -3*x +1"

print("\n\n\n*******Root finding by Newton Raphson Method*******\n")
  

  
def NewtonRaphson(s,y): 
    tol=1e-6
    Nmax=2000

    s1=sym.sympify(s)
    x=sym.symbols('x')
    f=sym.lambdify(x,s1,"numpy")



    f_diff=sym.diff(s1) 
    diff_f=sym.lambdify(x,f_diff,"numpy") 

    ratio = f(y)/diff_f(y) 
    iter=0
    while (abs(ratio) >= tol and iter <=Nmax):  
        ratio = f(y)/diff_f(y) 
        y = y - ratio 
        iter+=1
      
    print(f"The value of root is : {y} \n it took {iter} number of iterations to find the root with the accuracy of {f(y)}") 


#alpha=input("enter the guess for Newton Raphson method \n")
#alpha=float(alpha)
#NewtonRaphson(s,alpha)
