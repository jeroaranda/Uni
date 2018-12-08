function [x,iter] = mlbfgs(fname,x0,m)
% Método de búsqueda de línea donde la  dirección de descenso se calcula
% por medio de la actualización BFGS con memoria limitada.
%
% Análisis Aplicado
% ITAM
% 4 de diciembre de 2018

%In
% fname .- cadena con el nombre de la funcion.
% x0 .- vector columna de dimension n que es el punto inicial del 
%       proceso de minimización
% Out
% x.- aproximacion a un minimo local de fname.
%iter .- iteraciones que procesa el método.
%----------------------------------------------------------------------
tol = 1.e-06;                    % tolerancia a la norma del gradiente
x   = x0;                        % punto inicial
n   = length(x);                 % dimensión del problema
gx  = gradiente(fname,x);        % gradiente en el punto
%------------------------------------------------------------------------
maxiter = 100;                   % número máximo de iteraciones          
iter    = 0;                     % contador de las itercaiones
s=zeros(n,m);                   %Matriz  de vectores que almacenará si=xi+1-xi;
y=zeros(n,m);                      %Matriz de  vectores que almacenarán las diferencias de los gradientes.
gama=zeros(n,1);                   %Vector que almacena gamas, definidas por 1/st*y.
flag=true;                  %bandera para evitar que el método se indetermine.
%El método iterativo funciona mientras las condiciones de tolerancia se
%satisfazcan.


while ( norm(gx) > tol  && iter < maxiter&flag)
    
    %La dirección de descenso se calcula con el método de memoria limitada
    %BFGS en dónde se actualiza la matriz B cada m iteraciones. Para esto
    %se toma en cuenta las iteraciones módulo m. Nuestro primer vector q es
    %el gradiente e iterativamente se va calculando los nuevos qs y
    %almacenando las alphas.
    n1=mod(iter,m);
    alpha=zeros(n1);
    q=gx;
    for i=1:n1
        alpha(n1-(i-1))=gama(n1-(i-1))*s(:,n1-(i-1))'*q;
        q=q-alpha(n1-(i-1))*y(:,n1-(i-1));
    end
   %En esta parte se obtiene la primera r multiplicando nuestra B por q, en donde si las iteraciones
   %son menores que m se utiliza la identidad.  De otra forma se utiliza la
   %fórmula en los lineamientos del proyecto.
    if(iter<m)
        r=eye(n)*q;
    else
        r=gama(m)/(y(:,m)'*y(:,m))*eye(n)*q;
    end
    %Aquí se calcula iterativamente las r, ocupando las betas y las alphas
    %calculadas anteriormente.
    for i=1:n1
        betai=gama(i)*y(:,i)'*r;
        r=r+(alpha(i)-betai)*s(:,i);
    end
    %P es nuestra dirección de descenso, si esta es positiva utilizamos la
    %dirección de máximo descenso.
    p=-r;
    pend = gx'*p;
    if(pend<0)
        p=-gx;
    end
    %Posteriormente se hace un método de búsqueda de linea en donde se
    %interpola una función cuadrática en los puntos (x,fx) y (xp,fxp) y la
    %pendiente en x y se encuentra su mínimo. La interpolación se obtiene
    %ajustando un polinomio Q(x)=ax^2+bx+c, donde Q(0)=f(xk),
    %Q(norm(p))=f(xk+p) y Q'(0)=gradiente(f(xk))*p. Su mínimo se encuentra
    %en -b/2a.
    fx = feval(fname,x);
    fxp=feval(fname,x+p);
    long=norm(p);
    alpha=-pend*long/(2*(fxp-fx-pend*long));
    %Aquí se obtiene nuestro nuevo punto en la dirección de descenso y
    %escalado por alpha, nuestro minimo en la interpolación cuadrática.
    xp=x+alpha*p;
 
   
   %Aquí se actualizan nuestras matrices de vectores y escalares además de
   %actualizar el nuevo punto y el nuevo gradiente. Si gama por alguna
   %razón se indetermina el método se detiene.
    gxp = gradiente(fname,xp);
    s(:,n1+1) = xp - x ;
    y(:,n1+1) = gxp - gx;
    a=(s(:,n1+1)'*y(:,n1+1));
    if(a==0)
        flag=false;
    end
    gama(n1+1) = 1/a;
    disp(fx);
    x=xp;
    gx = gxp;
    iter = iter+1;
end
