% implementación del algoritmo mlbfgs para para la función DIXMAAN, con
% dimensión n. Ingresar vector incial múltiplo de tres. 
n=1500; 
m=17; 
if( mod(n, 3)==0 )
% Punto sugerido por Andrei.  
x0=2*ones(n,1);
t = cputime;
[sol]=mlbfgs('dixmaan', x0,m);
e = cputime-t;
fopt=feval('dixmaan',sol);
gfin=gradiente('dixmaan',sol);
X=sprintf('CPUTime: %d ; norma del gradiente: %d ; valor óptimo de la función: %d',t, norm(gfin),fopt );
disp(X)
else
    disp('Dimensiones incorrectas. Ingrese vector x0 de dimensión múltiplo positivo de 3')
end