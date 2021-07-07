function A = laplace(Nx)
    %Matriz diagonal:
        B = ones(Nx);
        B = tril(B,+1);
        B = triu(B,-1);
        B = ((B - eye(Nx)));
        B = (-2)*eye(Nx) + B;
     %Matriz identidad:
        I = eye(Nx);
        
%====================%Matriz Laplace final%================================

A = kron(B,I) + kron(I,B);

%====================%Condiciones de Neumann%=================================
%LIBRO PAG 589 = solo ultima fila condiciones de neumann
count = 0;

A(1 : Nx, Nx+1 : 2*Nx) = 2*eye(Nx);
A((Nx*Nx) - Nx + 1 : (Nx*Nx), (Nx*Nx) -2*Nx +1 : (Nx*Nx) - Nx) = 2*eye(Nx);
A(2, 1) = 2;
A(Nx - 1, Nx) = 2; 
A(Nx*Nx - Nx + 2, Nx*Nx - Nx + 1) = 2;
A(Nx*Nx - 1, Nx*Nx) = 2;

    for i = 2:Nx-3
    A(i*Nx + 1, i*Nx + 2) = 2;
    A((i+1)*Nx,(i+1)*Nx - 1) = 2;  
    end
    
%=========================%Condiciones de Neumann%========================================= 
 A = sparse(A);
end