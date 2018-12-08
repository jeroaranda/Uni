function [fx] =  freuroth(x)
% Extensión de la función de Freudenstein & Roth 
% con dimensión n par positivo. 
% In.
% x.- vector columna de dimensión n.
% Out
%fx.- número real.
%
% Análisis Aplicado
% ITAM
%  4 de diciembre de 2018
  
n=length(x);
fx=0; 
   
for i= 1:(n/2)
    fx=fx+(-13+x(2*i-1)+( (5-x(2*i))*x(2*i)-2)*x(2*i))^2+(-29+x(2*i-i)+( (x(2*i)+1)*x(2*i)-14)*x(2*i)   )^2;
end

end