%Actividad 1

%Analisis de sensibilidad
A = [6,4 ; 8,4 ; 3,3];
b = [40;40;20];
c = [300 ; 200 ];
[x0,zo,ban,iter,sensinfo] = mSimplexMax(A,b,c)

cambiob = [45;45;25];
[x1,z1,ban1,iter1,sensinfo1] = mSimplexMax(A,cambiob,c)