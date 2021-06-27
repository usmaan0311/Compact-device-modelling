from scipy.misc import derivative as diff

def Newton(f,ig):
	Nmax=2000000
	tol=1e-3

	print(f"{f}")
	

	if ( f(ig) == 0 ):
		print(f" {ig} is root of the function\n")
		return ig, 0

	else: 

		iter = 0
		y = ig
		error = 1.0
		h=1e-9

		while( (iter < Nmax) and (error > tol) ):
			fp = diff(f,y,h)
			
			y1 = y - (f(y)/fp)

			error = abs( f(y1) - f(y) )#/abs( f(y1) )
			y = y1
			iter+=1

		print(f"value of root is: {y}, with accuracy {error}, no of iterations are: {iter}, value of function is: {f(y)}\n ")
		return y, iter
