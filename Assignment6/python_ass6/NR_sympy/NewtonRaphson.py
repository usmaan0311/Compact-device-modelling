import sympy as sym

def Newton(f,ig):
	Nmax=2000000
	tol=1e-10

	print(f"{f}")
	s1=sym.sympify(f)
	print(f"{s1}")
	x=sym.symbols('x')
	f=sym.lambdify(x,s1,modules=['numpy'])

	if ( f(ig) == 0 ):
		print(f" {ig} is root of the function {f} \n")

	else:

		f_diff=sym.diff(s1) 
		diff_f=sym.lambdify(x,f_diff,"numpy",dtype=object) 

		iter = 0
		y = ig
		error = 1.0
		h=1e-6

		while( (iter < Nmax) and (error > tol) ):
			fp = diff_f(y)
			
			y1 = y - (f(y)/fp)

			error = abs( f(y1) - f(y) )/abs( f(y1) )
			y = y1
			iter+=1

		print(f"value of root is: {y}, with accuracy {error}, no of iterations are: {iter}, value of function is: {f(y)}\n ")
		return y, iter
