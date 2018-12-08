function [fx] = dixmaan(x)
% funciones de DIXMAANA - DIXMAANL,se vale de 8 par�metros (alpha a delta)
% k1 a k4, requiere que el vector de entrada sea de dimensi�n positiva, n, 
% n m�ltiplo de 3.4
% La implementaci�n de par�metros es la indexada con un E en el art�culo
% An Unconstrained Optimization Test Functions Collection de Andrei, N.
% (2008), p.155
% In.
% x.- vector columna de dimensi�n n.
% Out
% fx.- n�mero real.
%
% An�lisis Aplicado
% ITAM
%  4 de diciembre de 2018

alpha=1;
beta=0; 
gamma=0.125; 
delta=0.125; 
k1=1; 
k2=0; 
k3=0;
k4=1;  
n=length(x);
m=n/3; 


a=0; 
for i= 1:n
    a=a+x(i)^2*(i/n)^k1;
end
a=alpha*a; 


b=0; 
for i= 1:(n-1)
    b=b+x(i)^2*(x(i+1)+x(i+1)^2)*(i/n)^k2;
end 
b=beta*b; 

c=0; 
for i= 1:(2*m)
    c=c+x(i)^2*x(i+m)^4*(i/n)^k3;
end
c=gamma*c; 

d=0; 
for i= 1:m
    d=d+x(i)*x(i+2*m)^4*(i/n)^k4;
end
d=delta*d;

fx=1+a+b+c+d; 
    

end