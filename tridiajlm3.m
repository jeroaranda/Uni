% implementación del algoritmo mlbfgs para para la función TRIDIA, con
% dimensión n. Ingresar vector incial múltiplo de tres. 
n=1000; 
m=3; 
% Punto sugerido por Andrei.  
x0=ones(n,1);
t = cputime;
[sol]=mlbfgs('tridia', x0,m);
e = cputime-t;
fopt=feval('tridia',sol);
gfin=gradiente('tridia',sol);
X=sprintf('CPUTime: %d ; norma del gradiente: %d ; valor óptimo de la función: %d',t, norm(gfin),fopt );
disp(X)

