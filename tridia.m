function [fx] =  tridia(x)
% Funci�n tridia con entrada,x, de dimensi�n 2 o superior.
% In.
% x.- vector columna de dimensi�n n>=2.
% Out
%fx.- n�mero real.
%
% An�lisis Aplicado
% ITAM
%  4 de diciembre de 2018
  
n=length(x);
alpha=2; 
beta=1;
gamma=1; 
delta=1; 


fx=gamma*(delta*x(1)-1 )^2; 

   
for i= 2:n
    fx=fx+i*(alpha*x(i) -beta*x(i-1) )^2;
end

    

end