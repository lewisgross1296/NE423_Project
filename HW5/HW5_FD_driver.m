%% NE 423 HW5 ~Computational Part~
% Aysia Demby
% Michael Gerard
% Lewis Gross

clear; clc;

r_f = 0.005 ; %[m]
d_g = 0.00005; %[m]
d_c = 0.0006; %[m]
r_c = r_f + d_g + d_c; %[m]
r_g = r_f + d_g; %[m]
k_f = 3 ; %[W/m-K]
k_c = 17; %[W/m-K]
k_g = .25; %[W/m-K]
% T_S = 750; % [K]
T_cool = 580; % [K]
h_cool = 25000 ; % [W/m^2*K]
LHR = 20000 ; %[W/m]
Q = LHR/(pi*r_f^2) ; % LHR = Q*pi*r_f^2 [W/m^3]

%Number of Grid Points in Each Region
Nfs = [5, 10 , 20, 40 , 50];
Ngs = [5, 10 , 20, 30 , 40];
Ncs = [5, 10 , 15, 20 , 30];
Ns = Nfs + Ngs + Ncs - 2 ;

%step sizes between grid points
h_f = r_f ./ (Nfs - 1) ;
h_g = d_g ./ (Ngs - 1) ;
h_c = d_c ./ (Ncs - 1) ;
errors = zeros(1, length(Ns));
GS_iters = zeros(1, length(Ns));

for k = 1:length(Ns) %currently ends line 97, change? %channging to end at
    %line 56, and then starting another for loop? Hmmm....
    Nf = Nfs(k) ; hf = h_f(k) ;
    Ng = Ngs(k) ; hg = h_g(k) ;
    Nc = Ncs(k) ; hc = h_c(k) ;
    N = Ns(k) ;
    r = unique([(0:hf:r_f),(r_f:hg:r_g), (r_g:hc:r_c)]') ;
    fuel_grid = linspace(0,r_f,Nf) ;
    gap_grid = linspace(r_f,r_g,Ng) ;
    clad_grid = linspace(r_g,r_c,Nc) ;
    [Tvec, a_grid] = fuel_temp_analytical(LHR,r_f,d_g,d_c,k_f,k_g,k_c,h_cool, ...
        T_cool,fuel_grid, gap_grid, clad_grid);
    
    %initializing bvec
    bvec = zeros(N, 1) ;
    % RHS of FD conduction equation, AT=b
    % condition at r=0
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
    count = count + 1 ;fd
    
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
    bvec(N,1) = LHR/(2*pi*h_cool*r_c) + T_cool ;
    
    A = sparse(ivec,jvec,avec);
    tol = 1e-6;
    [T, counts] = GaussSeidel(A,bvec,tol);
    GS_iters(k)=counts ;
    %comparing analytical to numerical
    errors(k) = norm(T - Tvec,Inf) ;
    figure(k); plot(r,T,'bo',r,Tvec,'r')
    xlabel('Radius in Pellet [m]')
    ylabel('Temperature K')
    title(['Comparing Numerical and Analytical Solution for ',num2str(N),' grid points'])
    legend('numerical','analytical')
end

% average scheme h for error comparison
average_h = (h_f.*Nfs + h_g.*Ngs + h_c.*Ncs)./Ns;

% % number of iterations vs N grid points
figure(k+1);plot(Ns,GS_iters,'r-o')
xlabel('Number of Grid Points')
ylabel('Number of Iterations in GS Algorithm)')
title('Comparing GS Iterations Required vs Number of Grid Points')

% % showing quadratic convergence for central FD between grid spacing and
% % error
figure(k+2);plot(log(average_h),log(errors),'k-o')
xlabel('log(h)')
ylabel('log(errors(h))')
title('Log-Log Plot for Grid Spacing vs Error')

% Linear Polyfits to try and see slope of previous plot, slope =2 means
% quadratic convergence is true
coefs_h_ave = polyfit(log(average_h),log(errors),1);
coefs_h_f = polyfit(log(h_f),log(errors),1);
coefs_h_c = polyfit(log(h_c),log(errors),1);
coefs_h_g = polyfit(log(h_g),log(errors),1);
