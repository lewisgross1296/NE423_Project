function [T, r] = radial_solver(Q,TCO,r_f,r_g,r_c,k_f,k_g,k_c,Nf,Ng,Nc)
    % This function computes the radial temperature profile for the given set
    % of parameters in the function call
    % The output is the tempreature at each radial node and the radial grid
    % Q,TCO,r_f,r_g,r_c,k_f,k_g,k_c,Nf,Ng,Nc TODO
    
    hf = r_f/(Nf-1) ;
    hg = (r_g-r_f)/(Ng-1) ;
    hc = (r_c-r_g)/(Nc-1) ;
    N = Nf+Ng+Nc-2 ;
    LHR = Q*pi*r_f^2;
    %radial mesh
    r = unique([(0:hf:r_f),(r_f:hg:r_g), (r_g:hc:r_c)]') ;
    %initializing bvec
    bvec = zeros(N, 1) ;
    % RHS of FD conduction equation, AT=b
    % condition at r=0, Neumann BC
    bvec(1,1) = - hf^2/(2*k_f)*Q ;
    
    count = 1; % counter for sparse matrix formulation
    % Neumann Boundary, r=0
    ivec(count) = 1 ; jvec(count) = 1; avec(count) = - 1 ;
    count = count + 1 ;
    ivec(count) = 1 ; jvec(count) = 2 ; avec(count) =  1 ;
    count = count + 1 ;
    
    % Loop for internal nodes in the fuel
    for j = 2 : Nf - 1
        % j-1 (left of the diagonal)
        ivec(count) = j ; jvec(count) = j-1 ; avec(count) = r(j)*k_f/(hf^2) - k_f/(2*hf) ;
        count = count + 1 ;
        %j (diagonal)
        ivec(count) = j ; jvec(count) = j ; avec(count) = -2*r(j)*k_f/(hf^2) ;
        count = count + 1 ;
        %j+1 (right of the diagonal)
        ivec(count) = j ; jvec(count) = j+1 ; avec(count) = r(j)*k_f/(hf^2) + k_f/(2*hf) ;
        count = count + 1 ;
    end
    % heat generation in the fuel
    bvec(2: Nf-1,1) = - Q * r(2:Nf-1) ;
    
    %Continuity of Heat Flux Boundary Condition (fuel to gap)
    ivec(count) = Nf ; jvec(count) = Nf  ; avec(count) =  -1;
    count = count + 1 ;
    % right of diagonal
    ivec(count) = Nf  ; jvec(count) = Nf + 1; avec(count) = 1  ;
    count = count + 1 ;
    % heat flux enetering boundary
    bvec(Nf) = - LHR* hg / (2*pi*r_f*k_g) ;
    
    %Gap Region
    %Internal Nodes for gap region
    %Fuel gap boundary is at Nf, gap clad boundary is at Nf + Ng - 1
    for j = Nf + 1 : Nf + Ng - 2
        % j-1 (left of the diagonal)
        ivec(count) = j ; jvec(count) = j-1 ; avec(count) = r(j)*k_g/(hg^2) - k_g/(2*hg) ;
        count = count + 1 ;
        %j (diagonal)
        ivec(count) = j ; jvec(count) = j ; avec(count) = -2*r(j)*k_g/(hg^2) ;
        count = count + 1 ;
        %j+1 (right of the diagonal)
        ivec(count) = j ; jvec(count) = j+1 ; avec(count) = r(j)*k_g/(hg^2) + k_g/(2*hg) ;
        count = count + 1 ;
    end
    
    %Continuity of Heat Flux Boundary Condition (gap to clad)
    ivec(count) = Nf + Ng -1 ; jvec(count) = Nf + Ng -1  ; avec(count) =  -1;
    count = count + 1 ;
    % right of diagonal
    ivec(count) = Nf + Ng -1  ; jvec(count) = Nf + Ng; avec(count) = 1  ;
    count = count + 1 ;
    % heat flux enetering boundary
    bvec(Nf + Ng -1) = - LHR* hc / (2*pi*r_g*k_c) ;
    
    %Cladding Region
    %Internal Nodes for cladding region
    %gap clad boundary is at Nf + Ng -1 , clad coolant boundary is at Ns=
    %Nf + Ng + Nc - 2
    for j = Nf + Ng : Nf + Ng + Nc - 3
        % j-1 (left of the diagonal)
        ivec(count) = j ; jvec(count) = j-1 ; avec(count) = r(j)*k_c/(hc^2) - k_c/(2*hc) ;
        count = count + 1 ;
        %j (diagonal)
        ivec(count) = j ; jvec(count) = j ; avec(count) = -2*r(j)*k_c/(hc^2) ;
        count = count + 1 ;
        %j+1 (right of the diagonal)
        ivec(count) = j ; jvec(count) = j+1 ; avec(count) = r(j)*k_c/(hc^2) + k_c/(2*hc) ;
        count = count + 1 ;
    end
    
    % cladding to coolant boundary
    % Dirichlet Boundary
    ivec(count) = N ; jvec(count) = N ; avec(count) = 1;
    bvec(N,1) = TCO;
    
    A = sparse(ivec,jvec,avec);
    tol = 1e-6;
    [T, ~] = GaussSeidel(A,bvec,tol);
end