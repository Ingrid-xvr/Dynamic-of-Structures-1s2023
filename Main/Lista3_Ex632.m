clc
close all
clear all
% Definindo a matriz de massa M e matriz de rigidez K
M = [4,1,0,0;1,4,1,0;0,1,4,1;0,0,1,2];
K = [2,-1,0,0;-1,2,-1,0;0,-1,2,-1;0,0,-1,1];

% Initialization of variables
m = 0;
k = 1;

Phi1 = eye(size(M));
lambda = zeros(size(M, 1), 1); %n�o tenho certeza se � isso
lambda_old = lambda; %n�o tenho certeza se � isso

% Select threshold
epsilon = 10^(-8);

sizeM = size(M);
% Iteration until convergence
while true
    % Increment sweeping counter
    m = m + 1;
    
    % Sweeping
    for i = 1:sizeM-1 %n�o tenho certeza se � isso
        for j = i+1:sizeM %n�o tenho certeza se � isso
            % Test for rotation
            epsilon_r = 10^(-2*m);
            
            if sqrt((M(i,j))^2/((M(i,i))*(M(j,j)))) <= epsilon_r && sqrt((K(i,j))^2/(K(i,i)*(K(j,j)))) <= epsilon_r
                continue; % Sem rota��o, retorne ao loop externo
            end
            i;
            j;

            % Determination of transformation matrix Pk
            a0 = K(i,i)*M(i,j) - M(i,i)*K(i,j);
            a1 = K(j,j)*M(i,j) - M(j,j)*K(i,j);
            a2 = K(i,i)*M(j,j) - M(i,i)*K(j,j);
            
            aux = sign(a2);
            
            if a2 == 0
                aux = 1;
            end
            
            %a express�o sign(a2) � usada para determinar o sinal de a2 para fins de c�lculos adicionais
            gamma = a2/2 + aux*sqrt((a2/2)^2 + a0*a1); 
            
            if gamma == 0
                alpha = -K(i,j)/K(j,j);
                beta = 0;
            else
                alpha = -a0/gamma;
                beta = a1/gamma;
            end
            
            alpha;
            beta;
            
            % Rotation of matrix M1 and K
            Pk = eye(size(M));
            Pk(j,i) = alpha;
            Pk(i,j) = beta;
            
            Pk;
            M = Pk' * M * Pk;
            K = Pk' * K * Pk;
            
            % Update eigenvectors
            Phi1 = Phi1 * Pk;
            
            % Increment rotation counter
            k = k + 1;
        end
    end
    
    % Update eigenvalues
    for i = 1:sizeM %n�o tenho certeza
        lambda(i) = K(i,i) / M(i,i);
    end
    
    % Test for convergence
    if max(abs(lambda - lambda_old)./lambda) <= epsilon %n�o tenho certeza do max
        break; % Converged, exit the loop
    end
    
    lambda_old = lambda;
end

% Final evaluation
num_sweeps = m;
num_rotations = k - 1;

% Association of the approximate eigenvalues
lambda_approx = lambda;

% Association of the approximate eigenvectors
mass = eye(sizeM);

for i = 1 : sizeM
    mass(i,i) = 1/sqrt(M(i,i));
end

Phi_approx = Phi1*mass 

% Display results
disp("Number of sweeps: " + num_sweeps);
disp("Number of rotations: " + num_rotations);
disp("Eigenvalues:");
disp(lambda_approx);
disp("Eigenvectors:");
disp(Phi_approx);


