% implementaci�n del algoritmo mlbfgs para para la funci�n DIXMAAN, con
% dimensi�n n. Ingresar vector inicial par. 
n=1000; 
m=29; 
%punto inicial sugerido por Andrei: ( [0.5 2 0.5 2,...,0.5 2] )
if( mod(n, 2)==0 )
x0=0.5*ones(n,1);
for i=1:n
    if( mod(i,2)==0 )
        x0(i)=x0(i)-2.5; 
    end
end
t = cputime;
[sol]=mlbfgs('freuroth', x0,m);
e = cputime-t;

fopt=feval('freuroth',sol);
gfin=gradiente('freuroth',sol);
X=sprintf('CPUTime: %d ; norma del gradiente: %d ; valor �ptimo de la funci�n: %d',t, norm(gfin),fopt );
disp(X)
else
    disp('Dimensiones incorrectas. Ingrese vector x0 de dimensi�n m�ltiplo positivo de 3')
end