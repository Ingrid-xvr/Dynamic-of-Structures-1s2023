//Exercicio 22.2 - Verificacao no LAS

//Carga externa
Dbl P = Zero(1501,2)

// Matriz de Massa - concentrada nos pavimentos (20000 kg) 
Dbl M = [20000	0		
		0	20000	]

//Nao ha amortecimento
Dbl C = [0	0
	   0	0]

//Matriz de rigidez condensada
Dbl K = [36000000	     -18000000
	   -18000000	18000000]

//Descolocamento inicial
Dbl U0 = [0.02
	    0.02]

//Velocidade inicial
Dbl V0 = [0
	    0]

//Resposta da estrutura 
//U = Houbolt(P, M, C, K, U0, V0, delta_t) 
Dbl U = Houbolt(P, M, C, K, U0, V0, 0.01) 
Dbl Umax = Col_Max(Abs(U))

