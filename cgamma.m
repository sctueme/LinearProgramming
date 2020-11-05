function [ Gamma ] = cgamma(A_B, A_N, B , r_n)

%Encontramos intervalos maximales para constantes gamma tales que la
%solucion siga siendo optima al modificar las entradas del vector 
%de costos, c.
%Gamma es una matriz de tama?o nx2, y A_B, A_N, B ,r_n estan definidos 
%como en el metodo simplex revisado. 

H = (A_B) \ (A_N);
[m,n] = size(H);


%La matriz ii tiene en la j-esima fila, los posibles valores para el limite
%inferior de la j-esima gamma, id los tiene para el limite superior
ii = [];
id = [];

% a es un vector que contendra el limite inferior de la j-esima gamma, 
% b el limite superior
a = zeros(n,1);
b = zeros(n,1);
 
%Para cada una de las n j's tenemos dos opciones: 

for  j = 1 : n 
    
     %1.- j es un indice de la base: tenemos 3 posibles casos para Hjk

     if any(B == j) 
          
         aux = (B==j);
         for k = 1 : n
            %si Hjk es negativo, necesitamos cota inferior para gamma
            if H(aux,k) < 0
                ii(j,k) = (-r_n(k))/H(aux,k);
                id(j,k) = inf;
                
            %si Hjk es positivo, necesitamos cota superior para gamma
            else 
                if H(aux,k) > 0
                    ii(j,k) = -inf;
                    id(j,k) = (-r_n(k))/H(aux,k);
                    
                %si Hjk es 0, entonces gamma es libre
                else 
                    ii(j,k) = -inf;
                    id(j,k) = inf;
                end

             end
          a(j) = max(ii(j,:));
          b(j) = min(id(j,:));     
         end


      %si j no esta en la base, entonces la cota inferior de gamma es rk y no hay
      %cota superior
      else 
            for k = 1:n
            ii(j,k) = r_n(k);
            id(j,k) = inf;
            end 
        
      end

end 


%Como para este problema especifico necesitamos resolver maximizacion, 
%multiplicamos anteriormente el vector c por -1 y por ende tenemos que
%cambiar tambien los intervalos
Gamma = [-b,-a];





end

