using Calculus
#NewtonRaphson solver

function NewtonRaphson(f::Function,ig)
  Nmax=200000000
  tol=1e-12
  
 
  if (f(ig)==0)
    println(ig,"is the root of the function",f)

  else 

    iter=0
    y=ig
    error=1
    while( iter<Nmax && error>tol)
      fp=derivative(f,y)
    
      y1 = y - (f(y)/fp)

      error=abs(f(y1)-f(y))#/abs(f(y1))
      y=y1
      iter+=1

     end  
    println("value of root is: ",y," with accuracy of: ",error," no. of iterations it took are: ",iter," function value is: ",f(y));
    return [y, iter]
   end
end







  
  
