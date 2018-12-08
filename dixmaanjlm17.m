% implementaci�n del algoritmo mlbfgs para para la funci�n DIXMAAN, con
% dimensi�n n. Ingresar vector incial m�ltiplo de tres. 
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
X=sprintf('CPUTime: %d ; norma del gradiente: %d ; valor �ptimo de la funci�n: %d',t, norm(gfin),fopt );
disp(X)
else
    disp('Dimensiones incorrectas. Ingrese vector x0 de dimensi�n m�ltiplo positivo de 3')
end