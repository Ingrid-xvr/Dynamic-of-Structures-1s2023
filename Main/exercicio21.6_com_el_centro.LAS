//Exercicio 21.6 com El Centro

// Nos
Int XY = Zero(41,2)
Int aux = 0
For (Int i = 1 ; i <= 41 ; i++)
XY(i,2) = aux
aux = aux + 5
Next

// Conectividade dos elementos quanto aos nos
Int Elementos = Zero(40,2)
For (Int i = 1 ; i <= 40 ; i++)
Elementos(i,1) = i
Elementos(i,2) = i+1
Next
Int Beam_EL = Elementos

// DOF (x,y,theta)
Int DOF = Zero(41,3)
For (Int i = 2 ; i <= 41 ; i++)
DOF(i,1) = i-1
Next
Int EQ = DOF

// Locacao dos elementos
Int Beam_LM = Beam_Make_LM(Beam_EL,EQ) 

// Rigidez lateral por pavimento (k=12EI/l^3 em [N/m])
Dbl k = 10^8

// Comprimento do elemento (length)
Dbl length = 5// [m]

// Modulo de elasticidaede 
Dbl E = (k*length^3)/12// [N/m2] Inercia unitaria

// Propriedade dos materiais
Dbl PC = [1.0 0.0 1.0 E 0.0 0.0] // Coluna

// Matriz de rigidez
Dbl Ki = Zero(6,6)// Matriz de rigidez elemento analisado
Dbl K = Zero(40) // Matriz de rigidez global
For (Int i = 1 ; i <= 40 ; i++)
Ki = Beam_Make_K(i,XY,Beam_EL,PC,"Global")
Assemble(K,Ki,Beam_LM,i)
Next

// Massa por pavimento
Dbl m = 10^5 //[kg]
Dbl M = Zero(4)

// Matriz de massa
For (Int i = 1 ; i <= 39 ; i++)
M(i,i) = m
Next
M(40,40) = 7*m

// Matriz de amortecimento - nao ha
Dbl C = Zero(40,40)

// Analise modal
Dbl {Phi,W2} = Jacobi(M,K) // Normalizado pela matriz de massa
Dbl W = Sqrt(W2)
Dbl Periods = W@^(-1)@*2*3.141592654 

// Periodo critico
Dbl Tmenor = Col_Min(Periods)
Dbl Tcrit_dc = Tmenor/Pi //dif. centradas
Dbl Tcrit_nw = 0.551*Tmenor

// Fator de distribuicao da direcao da excitacao do terremoto quanto aos graus de liberdade(n)
Int r = Zero(40,1)
For (Int i = 1 ; i <= 40 ; i++)
r(i,1) = 1
Next

// Aplicacao terremoto
Input "ELCnorm.LAS"
Dbl ELCg = 0.349*9.81*ELCnorm 
Dbl Peff = -M*r*t(ELCg)

// Condicoes iniciais
Dbl u0 = Zero(40,1)
Dbl v0 = Zero(40,1)
Dbl a0 = Zero(40,1)

// Incrementos temporais utilizados (s)
Dbl delta_t1 = 0.005
Dbl delta_t2 = 0.01
Dbl delta_t3 = 0.02
Dbl delta_t4 = 0.04
Dbl delta_t5 = 0.08

// Metodos de integracao numerica:

// Diferencas centradas 
Dbl U_Central_dif_t1 = Central_Difference(t(Peff), M, C, K,u0,v0,a0,delta_t1) // t=0.005s
Dbl U_Central_dif_t2 = Central_Difference(t(Peff), M, C, K,u0,v0,a0,delta_t2) // t=0.01s
Dbl U_Central_dif_t3 = Central_Difference(t(Peff), M, C, K,u0,v0,a0,delta_t3) // t=0.02s
Dbl U_Central_dif_t4 = Central_Difference(t(Peff), M, C, K, u0,v0,a0,delta_t4) // t=0.1s
Dbl U_Central_dif_t5 = Central_Difference(t(Peff), M, C, K, u0,v0,a0,delta_t5) // t=1s
Dbl U_Central_dif_t1_40 = Col(U_Central_dif_t1,40)
Dbl U_Central_dif_t2_40 = Col(U_Central_dif_t2,40)
Dbl U_Central_dif_t3_40 = Col(U_Central_dif_t3,40)
Dbl U_Central_dif_t4_40 = Col(U_Central_dif_t4,40)
Dbl U_Central_dif_t5_40 = Col(U_Central_dif_t5,40)

// Newmark com aceleracao media
Dbl U_Newmark_avrg_t1 = Newmark_Average_Acceleration(t(Peff),M,C,K,u0,v0,delta_t1) // t=0.02s
Dbl U_Newmark_avrg_t2 = Newmark_Average_Acceleration(t(Peff),M,C,K,u0,v0,delta_t2) // t=0.02s
Dbl U_Newmark_avrg_t3 = Newmark_Average_Acceleration(t(Peff),M,C,K,u0,v0,delta_t3) // t=0.02s
Dbl U_Newmark_avrg_t4 = Newmark_Average_Acceleration(t(Peff),M,C,K,u0,v0,delta_t4) // t=0.03s
Dbl U_Newmark_avrg_t5 = Newmark_Average_Acceleration(t(Peff),M,C,K,u0,v0,delta_t5) // t=0.10s
Dbl U_Newmark_avrg_t1_40 = Col(U_Newmark_avrg_t1,40)
Dbl U_Newmark_avrg_t2_40 = Col(U_Newmark_avrg_t2,40)
Dbl U_Newmark_avrg_t3_40 = Col(U_Newmark_avrg_t3,40)
Dbl U_Newmark_avrg_t4_40 = Col(U_Newmark_avrg_t4,40)
Dbl U_Newmark_avrg_t5_40 = Col(U_Newmark_avrg_t5,40)

//Newmark com aceleracao linear
Dbl U_Newmark_lin_t1 = Newmark_Linear_Acceleration(t(Peff), M, C, K, u0,v0,delta_t1) // t=0.02s
Dbl U_Newmark_lin_t2 = Newmark_Linear_Acceleration(t(Peff), M, C, K, u0,v0,delta_t2) // t=0.02s
Dbl U_Newmark_lin_t3 = Newmark_Linear_Acceleration(t(Peff), M, C, K, u0,v0,delta_t3) // t=0.02s
Dbl U_Newmark_lin_t4 = Newmark_Linear_Acceleration(t(Peff), M, C, K, u0,v0,delta_t4) // t=0.03s
Dbl U_Newmark_lin_t5 = Newmark_Linear_Acceleration(t(Peff), M, C, K, u0,v0,delta_t5) // t=0.10s
Dbl U_Newmark_lin_t1_40 = Col(U_Newmark_lin_t1,40)
Dbl U_Newmark_lin_t2_40 = Col(U_Newmark_lin_t2,40)
Dbl U_Newmark_lin_t3_40 = Col(U_Newmark_lin_t3,40)
Dbl U_Newmark_lin_t4_40 = Col(U_Newmark_lin_t4,40)
Dbl U_Newmark_lin_t5_40 = Col(U_Newmark_lin_t5,40)


// Wilson Theta
Dbl U_WT_t1 = Wilson_Theta(t(Peff), M, C, K, u0,v0,delta_t1, 1.41) // t=0.02s
Dbl U_WT_t2 = Wilson_Theta(t(Peff), M, C, K, u0,v0,delta_t2, 1.41) // t=0.02s
Dbl U_WT_t3 = Wilson_Theta(t(Peff), M, C, K, u0,v0,delta_t3, 1.41) // t=0.02s
Dbl U_WT_t4 = Wilson_Theta(t(Peff), M, C, K, u0,v0,delta_t4, 1.41) // t=0.03s
Dbl U_WT_t5 = Wilson_Theta(t(Peff), M, C, K, u0,v0,delta_t5, 1.41) // t=0.10s
Dbl U_WT_t1_40 = Col(U_WT_t1,40)
Dbl U_WT_t2_40 = Col(U_WT_t2,40)
Dbl U_WT_t3_40 = Col(U_WT_t3,40)
Dbl U_WT_t4_40 = Col(U_WT_t4,40)
Dbl U_WT_t5_40 = Col(U_WT_t5,40)

// HHT-alpha
Dbl alpha = -1/3
Dbl U_HHT_t1 = HHT_Alpha(t(Peff), M, C, K, u0,v0,delta_t1, alpha) // t=0.02s
Dbl U_HHT_t2 = HHT_Alpha(t(Peff), M, C, K, u0,v0,delta_t2, alpha) // t=0.02s
Dbl U_HHT_t3 = HHT_Alpha(t(Peff), M, C, K, u0,v0,delta_t3, alpha) // t=0.02s
Dbl U_HHT_t4 = HHT_Alpha(t(Peff), M, C, K, u0,v0,delta_t4, alpha) // t=0.03s
Dbl U_HHT_t5 = HHT_Alpha(t(Peff), M, C, K, u0,v0,delta_t5, alpha) // t=0.10s
Dbl U_HHT_t1_40 = Col(U_HHT_t1,40)
Dbl U_HHT_t2_40 = Col(U_HHT_t2,40)
Dbl U_HHT_t3_40 = Col(U_HHT_t3,40)
Dbl U_HHT_t4_40 = Col(U_HHT_t4,40)
Dbl U_HHT_t5_40 = Col(U_HHT_t5,40)

// Houbolt
Dbl U_houbolt_t1 = Houbolt(t(Peff), M, C, K, u0,v0,delta_t1) // t=0.02s
Dbl U_houbolt_t2 = Houbolt(t(Peff), M, C, K, u0,v0,delta_t2) // t=0.02s
Dbl U_houbolt_t3 = Houbolt(t(Peff), M, C, K, u0,v0,delta_t3) // t=0.02s
Dbl U_houbolt_t4 = Houbolt(t(Peff), M, C, K, u0,v0,delta_t4) // t=0.03s
Dbl U_houbolt_t5 = Houbolt(t(Peff), M, C, K, u0,v0,delta_t5) // t=0.10s
Dbl U_houbolt_t1_40 = Col(U_houbolt_t1,40)
Dbl U_houbolt_t2_40 = Col(U_houbolt_t2,40)
Dbl U_houbolt_t3_40 = Col(U_houbolt_t3,40)
Dbl U_houbolt_t4_40 = Col(U_houbolt_t4,40)
Dbl U_houbolt_t5_40 = Col(U_houbolt_t5,40)

//Grau de liberdade 40 para cada periodo com os diferentes metodos
Dbl U40_aux1 = Combine(Combine(U_Central_dif_t1_40,U_Newmark_avrg_t1_40),Combine(U_Newmark_lin_t1_40,U_WT_t1_40))
Dbl U40_t1 = Combine(U40_aux1,Combine(U_HHT_t1_40,U_houbolt_t1_40))

Dbl U40_aux2 = Combine(Combine(U_Central_dif_t2_40,U_Newmark_avrg_t2_40),Combine(U_Newmark_lin_t2_40,U_WT_t2_40))
Dbl U40_t2 = Combine(U40_aux2,Combine(U_HHT_t2_40,U_houbolt_t2_40))

Dbl U40_aux3 = Combine(Combine(U_Central_dif_t3_40,U_Newmark_avrg_t3_40),Combine(U_Newmark_lin_t3_40,U_WT_t3_40))
Dbl U40_t3 = Combine(U40_aux3,Combine(U_HHT_t3_40,U_houbolt_t3_40))

Dbl U40_aux4 = Combine(Combine(U_Central_dif_t4_40,U_Newmark_avrg_t4_40),Combine(U_Newmark_lin_t4_40,U_WT_t4_40))
Dbl U40_t4 = Combine(U40_aux4,Combine(U_HHT_t4_40,U_houbolt_t4_40))

Dbl U40_aux5 = Combine(Combine(U_Central_dif_t5_40,U_Newmark_avrg_t5_40),Combine(U_Newmark_lin_t5_40,U_WT_t5_40))
Dbl U40_t5 = Combine(U40_aux5,Combine(U_HHT_t5_40,U_houbolt_t5_40))

//t4 que convergiram
Dbl U40_aux6 = Combine(Combine(U_Newmark_avrg_t4_40,U_Newmark_lin_t4_40),Combine(U_WT_t4_40,U_HHT_t4_40))
Dbl U40_t4converg = Combine(U40_aux6,U_houbolt_t4_40)

//t5 que convergiram
Dbl U40_t5converg = Combine(Combine(U_Newmark_avrg_t5_40,U_WT_t5_40),Combine(U_HHT_t5_40,U_houbolt_t5_40))



