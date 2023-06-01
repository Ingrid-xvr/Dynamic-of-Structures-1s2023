//Exercicio 21.4

//Matriz de rigidez condensada
Dbl K = [25335200	      -14381600	3344280	-427270
         -14381600	22084600	-13666400	2665920
           3344280	-13666400	20654300	-9821180
           -427270	2665920	-9821180	7517250]

// Matriz de Massa - concentrada nos pavimentos (20000 kg) 
Dbl M = [20000	0		0		0
		0	20000		0		0
		0	0		20000		0
		0	0		0		20000]

//Modal analysis
Dbl {Phi,W2} = Jacobi(M,K)
Dbl W = Sqrt(W2)
Dbl Periods = W@^(-1) @* 2*3.141592654

//Inverse Iteration Problem with Spectral Shift
Dbl Delta1 = 1938 //Delta é um fator de decalagem espectral, fica chutndo valor para achar o W1, W2...
Dbl e1 = 1
Dbl X = Init(4,1,1) //vetor unitario 4x1
Dbl Y = M*X
Int k = 0
Dbl U = 0
Dbl fact
Dbl fact1

While ( e1 > 1.0E-6 ) //precisao de 10^-6
	k = k+1
	Dbl x2 = Inv(K-Delta1*M)*Y
	Dbl y2 = M*x2
	Dbl X_t = Transpose(x2)
	Dbl Delta2 = ((X_t*Y)/(X_t*y2)) + Delta1
	fact = Sqrt(X_t*y2)
	fact1 = 1/fact
	U = fact1*y2
	e1 = Abs(Delta2-Delta1)@/Delta2
	Delta1 = Delta2
	Y = U

EndWhile

Dbl W_inverse_iteration = Sqrt(Delta1)
Dbl PHI_inverse_iteration = Inv(M)*U
End