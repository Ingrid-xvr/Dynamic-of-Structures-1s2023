//Exercicio 18.1

//Nos do portico
Int XY = [ 0	0    
	     3	0             
           0	3    
           3	3    
           0	6    
	     3	6   
	     0	9    
	     3	9    
	     0	12    
	     3	12]    
            
//Elementos do portico para serem unidos depois
Int Beam_EL = [1	3 //pilar esquerda 1 pav       
		   2	4 //pilar direita 1 pav
 		   3	4 //viga 1 pav
		   3	5 //pilar esquerda 2 pav
 		   4	6 //pilar direita 2 pav
		   5	6 //viga 2 pav
		   5	7 //pilar esquerda 3 pav
		   6	8 //pilar direita 3 pav
		   7	8 //viga 3 pav
		   7	9 //pilar esquerda 4 pav
		   8	10 //pilar direita 4 pav
		   9	10] //viga 4 pav  
               
//DOF's - translacional na horizontal (x,y,rotacao)
Int EQ = [ 0 0 0    //no 1
           0 0 0    //no 2          
           1 0 0    //no 3
           1 0 0    //no 4
           2 0 0    //no 5
           2 0 0    //no 6
           3 0 0    //no 7
           3 0 0    //no 8
           4 0 0    //no 9
           4 0 0]   //no 10
         			 
//Locacao das vigas na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL,EQ)

// Properties matrices [A As I E Nu Dn]
Dbl PC = [1.0 0.0 1.0 16875000 0.0 0.0]
Dbl PB = [1.0 0.0 1.0 1.0 0.0 0.0]

// Beams stiffness matrices
Dbl KC1 = Beam_Make_K(1,XY,Beam_EL,PC,"Global")
Dbl KC2 = Beam_Make_K(2,XY,Beam_EL,PC,"Global")
Dbl KB1 = Beam_Make_K(3,XY,Beam_EL,PB,"Global")
Dbl KC3 = Beam_Make_K(4,XY,Beam_EL,PC,"Global")
Dbl KC4 = Beam_Make_K(5,XY,Beam_EL,PC,"Global")
Dbl KB2 = Beam_Make_K(6,XY,Beam_EL,PB,"Global")
Dbl KC5 = Beam_Make_K(7,XY,Beam_EL,PC,"Global")
Dbl KC6 = Beam_Make_K(8,XY,Beam_EL,PC,"Global")
Dbl KB3 = Beam_Make_K(9,XY,Beam_EL,PB,"Global")
Dbl KC7 = Beam_Make_K(10,XY,Beam_EL,PC,"Global")
Dbl KC8 = Beam_Make_K(11,XY,Beam_EL,PC,"Global")
Dbl KB4= Beam_Make_K(12,XY,Beam_EL,PB,"Global")

// Matriz de Rigidez Global 4x4 (4 DOF's)
Dbl K = Zero(4)
Assemble(K,KC1,Beam_LM,1) //elemento 1
Assemble(K,KC2,Beam_LM,2) //elemento 2
Assemble(K,KB1,Beam_LM,3) //elemento 3
Assemble(K,KC3,Beam_LM,4) //elemento 4 
Assemble(K,KC4,Beam_LM,5) //elemento 5
Assemble(K,KB2,Beam_LM,6) //elemento 6 
Assemble(K,KC5,Beam_LM,7) //elemento 7
Assemble(K,KC6,Beam_LM,8) //elemento 8
Assemble(K,KB3,Beam_LM,9) //elemento 9
Assemble(K,KC7,Beam_LM,10) //elemento 10
Assemble(K,KC8,Beam_LM,11) //elemento 11 
Assemble(K,KB4,Beam_LM,12) //elemento 12

// Matriz de Massa - concentrada nos pavimentos (80000 kg) 
Dbl M = [40000	0		0		0
		0	40000		0		0
		0	0		40000		0
		0	0		0		40000]

// Analise modal
Dbl {Phi,W2} = Jacobi(M,K)
Dbl W = Sqrt(W2)
Dbl Periods= W@^(-1) @* 2*3.141592654

Dbl Mgen = t(Phi)*M*Phi //Massa modal (ou massa generalizada para normalizar Phi)
Dbl Kgen = t(Phi)*K*Phi //Rigidez modal (rigidez generalizada)

// Direcao da excitacao do terremoto 
Dbl r = [1;1;1;1] // Graus de liberdade 1,2,3,4 e excitacao - horizontal

// Fator de participacao Modal
Dbl L = t(Phi)*M*r

// Vetor aceleracao
Input "ELCnorm.LAS"
Dbl ELCg = 0.349*9.81*ELCnorm 

// Taxa de amortecimento critico para cada modo (Xi=5%)
Dbl Xi=[0.05
	  0.05
	  0.05
	  0.05]

// Velocidades e deslocamentos iniciais
Dbl U0 = Zero(4,1)
Dbl V0 = Zero(4,1)

// Piecewuise Linear - retorna os vetores de deslocamento e velocidade para os sistemas desacoplados
// usando o esquema de integração linear por partes - coordenadas generalizadas
Dbl Z_modal = Piecewise_Linear(ELCg, W, Xi, L, U0, V0, 0.02) //Delta_t=0.02 s (El Centro)

// Deslocamento total (coordenadas reais)
Dbl U_modal = t(Phi * t(Z_modal))

// Deslocamentos maximos - em modulo (analisa valores positivos e negativos)
Dbl Zmax = Col_Max(Z_modal)
Dbl Zmin = Col_Min(Z_modal)
Dbl Umax = Col_Max(U_modal)
Dbl Umin = Col_Min(U_modal)

// Forca elastica 
Dbl fs_modal = t(K * t(U_modal))

// Valores maximos
Dbl fsmax = Col_Max(fs_modal)
Dbl fsmin = Col_Min(fs_modal)
Dbl fs = Col_Max(Abs(fs_modal))

// Momento de tombamento (base da estrutura)
Dbl X = 3*[1; 2; 3; 4]
Dbl M0_modal = t(t(X)*t(fs_modal))

// Valores maximos (modulo)
Dbl M0max = Col_Max(M0_modal)
Dbl M0min = Col_Min(M0_modal)
Dbl M0 = Col_Max(Abs(M0_modal))
