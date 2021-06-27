include("NewtonRaphson.jl")


using SymPy

x=Sym("x")

η=1
I0=50e-9
R=1000
vdd=2
T=300
vt=0.026*T/(300);

f(x)=vdd - I0*(exp(x/(η*vt)) -1)*R - x 

y=NewtonRaphson(f,0);

I(x)=I0*(exp(x/(η*vt)) -1)
I(y)*1000

println("Current in the circuit is: ",I(y)*1000," mA")


