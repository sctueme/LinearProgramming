function [xo, zo, ban, iter, lamo] = mSimplexDual(A, b, c)
%   purpose: Versi??n del Simplex Dual (revisado)
%        min { c^T x }
%        sujeto a
%            Ax >= b,
%            x >= 0,
%            c >= 0
% In :
%     A ... m x n matrix
%     b ... column vector with as many rows as A
%     c ... column vector with as many columns as A
%
% Out: xo ... SFB ??ptima del problema
%      zo ... valor ??ptimo del problema
%      ban ... indica casos
%             -1 ... si el conjunto factible es vacio
%             0  ... si se encontro una soluci ??on  ??optima
%             1 ... si la funci ??on objectivo no es acotada.
%       iter ... es el numer ??o de iteraciones (cambios de variables basicas) quehizoelm ??etodo
%       lamo ... Soluci??n del problema dual

    [m,n] = size(A);
    A = deal([A -eye(m)]);
    c = deal([c zeros(1,m)]);
    n = n+m;
    B = [n-m+1:1:n];
    iteraciones=0;
    while iteraciones < 500000
      ix = 1;
      for(i=1:n)
        if(any(B(:)==i))
        else
          N(ix) = i;
          ix = ix+1;
        end
      end
      x_B = A(:,B)\b;
      if(min(x_B) >= 0)
          x = zeros(1,n-m);
          for(i=1:1:m)
            if(B(i)<=n-m)
              x(B(i)) = x_B(i);
            end
          end
          xo = x;
          zo = c(B)*x_B;
          ban = 0;
          sol_dual = A(:,B)'\c(B)';
          iter = iteraciones;
          lamo = sol_dual;
          return;
        end
        salida = -1;
        for(i=1:1:m)
          if(x_B(i)<0)
              salida = i;
              break;
          end
        end
        H = A(:,B)\A;
        %disp(c(B));
        %disp(H(:,N));
        r_N = -c(B)*H(:,N) + c(N);
        %disp(r_N)
        entrada = -1;
        for (i=1:1:n-m)
          if H(salida,N(i))<0
            if entrada == -1
               entrada = i;
            elseif r_N(i)/H(salida,N(i)) > r_N(entrada)/H(salida,N(entrada))
               entrada=i;
            end
          end
        end
        if entrada == -1
          xo = deal([]);
          zo = deal([]);
          ban = -1;
          iter = iteraciones;
          lamo = deal([]);
          return
        end
        B(salida)=N(entrada);
        iteraciones = iteraciones+1;
      end
end
