//Exercicio 18.2

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

// Matriz de Massa - concentrada nos pavimentos (40000 kg) 
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

Dbl Phi1 = Col(Phi,1)
Dbl Phi2 = Col(Phi,2)
Dbl Phi3 = Col(Phi,3)
Dbl Phi4 = Col(Phi,4)


// Pseudo displacement spectra - grafico
Dbl sd1 = 0.119
Dbl sd2 = 0.0178
Dbl sd3 = 0.0072
Dbl sd4 = 0.0050

// Displacement combinaison
Dbl U1M = Phi1*L(1,1)*sd1
Dbl U2M = Phi2*L(2,1)*sd2
Dbl U3M = Phi3*L(3,1)*sd3
Dbl U4M = Phi4*L(4,1)*sd4

Dbl UMax = SRSS(Combine(Combine(U1M,U2M),Combine(U3M,U4M)))

// Pseudo acceleration spectra - grafico
Dbl g=9.81
Dbl sa1 = 0.550*g
Dbl sa2 = 0.655*g
Dbl sa3 = 0.695*g
Dbl sa4 = 0.700*g

// Force combinaison
Dbl FS1M = M*Phi1*L(1,1)*sa1
Dbl FS2M = M*Phi2*L(2,1)*sa2
Dbl FS3M = M*Phi3*L(3,1)*sa3
Dbl FS4M = M*Phi4*L(4,1)*sa4

Dbl F=Combine(Combine(FS1M,FS2M),Combine(FS3M,FS4M))
Dbl FSM = SRSS(Combine(Combine(FS1M,FS2M),Combine(FS3M,FS4M)))

Dbl V11 = FS1M(1,1)+FS1M(2,1)+FS1M(3,1)+FS1M(4,1)
Dbl V21 = FS1M(2,1)+FS1M(3,1)+FS1M(4,1)
Dbl V31 = FS1M(3,1)+FS1M(4,1)
Dbl V41 = FS1M(4,1)
Dbl Vmax1=t Combine(Combine(V11,V21),Combine(V31,V41))

Dbl V12 = FS2M(1,1)+FS2M(2,1)+FS2M(3,1)+FS2M(4,1)
Dbl V22 = FS2M(2,1)+FS2M(3,1)+FS2M(4,1)
Dbl V32 = FS2M(3,1)+FS2M(4,1)
Dbl V42 = FS2M(4,1)
Dbl Vmax2=t Combine(Combine(V12,V22),Combine(V32,V42))

Dbl V13 = FS3M(1,1)+FS3M(2,1)+FS3M(3,1)+FS3M(4,1)
Dbl V23 = FS3M(2,1)+FS3M(3,1)+FS3M(4,1)
Dbl V33 = FS3M(3,1)+FS3M(4,1)
Dbl V43 = FS3M(4,1)
Dbl Vmax3=t Combine(Combine(V13,V23),Combine(V33,V43))

Dbl V14 = FS4M(1,1)+FS4M(2,1)+FS4M(3,1)+FS4M(4,1)
Dbl V24 = FS4M(2,1)+FS4M(3,1)+FS4M(4,1)
Dbl V34 = FS4M(3,1)+FS4M(4,1)
Dbl V44 = FS4M(4,1)
Dbl Vmax4 =t Combine(Combine(V14,V24),Combine(V34,V44))

Dbl Vmax = Combine(Combine(Vmax1,Vmax2),Combine(Vmax3,Vmax4))
Dbl Vmaxfinal = SRSS(Combine(Combine(Vmax1,Vmax2),Combine(Vmax3,Vmax4)))



