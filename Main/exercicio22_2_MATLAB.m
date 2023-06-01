clear all
%Dados de entrada
M = [20000,0;0,20000];
K = [36000000,-18000000;-18000000,18000000];

%Cálculos preliminares

%1. Cálculo das constantes de integração
dt = 0.01 %dado do exemplo 22.1

a0=2/dt^2;
a1=11/(6*dt);
a2=5/dt^2;
a3=3/dt;
a4=-4/dt^2;
a5=-3/(2*dt);
a6=1/dt^2;
a7=2/(3*dt);

%2.Determinação de u1 e u2 usando um procedimento especial. (Dado do exemplo 22.1)
u0 = [0.02; 0.02] %tempo zero
    
u1 = [0.0191; 0.02] %tempo 1

u2 = [0.016562; 0.019919] %tempo 2

%3. Determinação da matriz de rigidez efetiva
Kef = K + a0*M ;%não tem amortecimento 

%Integração passo a passo: para n = 2, 3,...,td/delta_t
%1. Tempo de incremento
td = 15 ;%adotado
nf = td/dt - 1;

un = zeros(2, td/dt+1);
un(:,1) = u0 %substituir todas as linhas da primeira coluna pelo vetor u0
un(:,2) = u1; %substituir todas as linhas da segunda coluna pelo vetor u1
un(:,3) = u2; %substituir todas as linhas da terceira coluna pelo vetor u2

t_inc = [0:dt:td]; %tempo analisado
p_ef = zeros(2, td/dt+1);
p_ef(:,2) = transpose(u1)*Kef;
p_ef(:,3) = transpose(u2)*Kef;

%2. Determinação do vetor de forças efetivas ˆpn+1 no tempo tn+1 e 3. Determinação do vetor de deslocamentos no tempo tn+1:

for jj = 2:nf ;
   p_ef(:,jj+2) = M*(a2*un(:,jj+1)+a4*un(:,jj)+a6*un(:,jj-1));
   un(:,jj+2) = transpose(p_ef(:,jj+2))*inv(Kef);
end

un_f = un'

figure(1)
 plot(t_inc, un(1,:))
 hold on
 plot(t_inc, un(2,:))
 xlabel('t [s]')
 ylabel('u [m]')

