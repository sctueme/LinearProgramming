function [x0, z0, ban, iter, sensinfo] = mSimplexMax(A, b, c)
% purpose: Versi?on del Simplex (revisado)
%   maximizar c^T x
%   sujeto a Ax <= b , x >= 0 , b >= 0
%
% In : A ... m x n matrix
%       b ... column vector with as many rows as A
%       c ... column vector with as many columns as A
%
% Out:  xo ... SFB ?optima del problema
%       zo ... valor ?optimo del problema
%       ban ... indica casos:
%           -1 ... si el conjunto factible es vacio
%            0 ... si se encontro una soluci?on ?optima
%            1 ... si la funci?on objectivo no es acotada.
%       iter ... es el numer?o de iteraciones (cambios de variables basicas)
%                que hizo el m?etodo
%
%       sensinfo ... Solo cuando ban = 0:
%                    sensinfo.lambda ... es la soluci?on dual
%                    sensinfo.gammas ... nx2 matrix con intervalos
%                    sensinfo.betas ... mx2 matrix con intervalos
%


    %Vemos el tama???o de las variables b???sicas y no b???sicas.
    [tam_B, tam_N] = size(A);
    N = 1:tam_N;
    B = (tam_N+1):(tam_N+tam_B);
    %Agregamos las variables de holgura
    A = [A eye(tam_B)];
    %Insertamos los ceros en los costos de las variables de holgura
    c = -[c' zeros(1,tam_B)];
    %Definimos las variables que faltan
    r_n = -c;
    h = b;
    iter = 0;
    %Vemos si algun valor del lado derecho es negativo. Si lo es, el
    %conjunto factible es vac???o y ya acabamos
    if any(h < 0)
        ban = -1;
    else
        ban = 0;
        %Ahora que ya sabemos que el problema es acotado, vemos si alg???n
        %costo relativo es mayor a cero
        while any(r_n > 0) && ban == 0
            %Utilizamos la Regla de Bland para definir las nuevas variables
            %b???sicas y no b???sicas
            entrada = find(r_n > 0, 1);
            N_e = N(entrada);
            col_bas = A(:, B)\A(:, N_e);
            %La variable asociada a la columna nueva_b entra como b???sica.
            %Preguntamos si la columna de la variable que entra es mayor a
            %cero, ya que sino el CF es no acotado y terminamos
            if any(col_bas > 0)
                iter = iter + 1;
                %Ahora veamos qu??? variable sale. No contamos las variables
                %de la columna que sean <= 0.
                %col_bas(col_bas<=0) = 0;
                cc = h./col_bas;
                cc(col_bas<=0) = inf;
                salida = find(min(cc)==cc);
                %Ahora que ya sabemos que variables entran y salen,
                %actualizamos los conjuntos B y N
                 B_s = B(salida);
                 B(salida) = N_e;
                 N(entrada) = B_s;
                %Obtenemos h
                h = A(:,B)\b;
                %Obtenemos los nuevos costos relativos
                r_n=((A(:,B)'\c(B)')'*A(:,N))-c(N);

            else
                %Como es no acotado, mandamos la bandera
                ban = 1;
            end
        end
    end
    %Ahora checamos si el problema es acotado, CF vac???o o con soluci???n
    %Si es no acotado o CF vac???o entonces devolvemos inf pues no hay sol.
    if ban == 1 || ban == -1
        x0 = inf;
        z0 = inf;
    else
        %Si hay soluci???n entonces obtenemos los valores de
        %x0 y z0
        x0 = zeros(tam_B+tam_N,1);
        x0(B) = h;
        z0 = -c*x0;
        x0 = x0(1:tam_N);


%Agregamos la informacion sobre sensibilidad.
        lambda = A(:,B)'\c(B)';
        sensinfo.lambda = -lambda';
        sensinfo.gammas = cgamma(A(:,B),A(:,N),B, r_n)
        sensinfo.betas = cbeta(A(:,B),h)
        gammas = cgamma(A(:,B),A(:,N),B, r_n);
        Betas = cbeta(A(:,B),h);
        

    end
    return
end