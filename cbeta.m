function [ Beta ] = cbeta(A_B, h)

%Aunque A_B es cuadrada (pues es invertible) utilizamos la siguiente linea 
%para obtener m explicitamente.
[m,l] = size(A_B);

%La matriz ii tiene en la i-esima fila, los posibles valores para el limite
%inferior de la i-esima beta, id los tiene para el limite superior
ii = zeros(m);
id = zeros(m);

% a es un vector que contendra el limite inferior de la i-esima beta, 
% b el limite superior
a = zeros(m,1);
b = zeros(m,1);

%Necesitaremos trabajar con la inversa de A_B;
%Matlab especifica que realizar la operacion * M\ * es mas eficiente que
%calcular Inv (inversa):
A = A_B\eye(m); 


for  k = 1 : m 
    
        for i = 1 : m
    
            %Tendremos 3 posibles casos para Aik: 
            %Lo siguiente se sigue de despejar analiticamente, como en
            %la funcion 'cgamma.m'
            if A(i,k) < 0
            ii(i,k) = -inf;
            id(i,k) = (-h(i))/A(i,k);
     
            else 
                if A(i,k) > 0
                ii(i,k) = (-h(i))/A(i,k);
                id(i,k) = inf;
  
                else 
                ii(i,k) = -inf;
                id(i,k) = inf;
                end
            end
        end

end 

%dado que cada beta necesita cumplir todas las cotas simultaneamente, 
%obtenemos el conjunto de la interseccion de posibles intervalos para
%beta:
for k = 1:m

    a(k) = max(ii(:,k));
    b(k) = min(id(:,k));    

end


Beta = [a,b];






end