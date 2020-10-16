%% NE 423 HW4 ~Computational Part~
% Aysia Demby
% Michael Gerard
% Lewis Gross

% clear; clc;
r_f = 0.005 ; %[m]
k_f = 3 ; %[W/m-K]
T_S = 750; % [K]
LHR = 20000 ; %[W/m]
Q = LHR/(pi*r_f^2) ; % LHR = Q*pi*r_f^2 [W/m^3]
Ns = [5 , 15, 45 , 90 ] ;
h = r_f ./ (Ns-1) ;
errors = zeros(1, length(Ns));
GS_iters = zeros(1, length(Ns));
for k = 1:length(Ns)
N = Ns(k);
r = [0:h(k):r_f]' ;
T_f = fuel_temp_analytical(r_f,k_f,LHR,T_S,r);
% RHS of FD conduction equation, AT=b
bvec = - Q * r ;
% condition at r=0
bvec(1) = - h(k)^2/(2*k_f)*Q ;
% condition at r=r_f
bvec(N) = T_S ;

count = 1; % counter for sparse matrix formulation
% Neumann Boundary, r=0
ivec(count) = 1 ; jvec(count) = 1; avec(count) = - 1 ;
count = count + 1 ;
ivec(count) = 1 ; jvec(count) = 2 ; avec(count) =  1 ;
count = count + 1 ;

% Loop for internal nodes
for j = 2 : N-1
% j-1 (left of the diagonal)
ivec(count) = j ; jvec(count) = j-1 ; avec(count) = r(j)*k_f/(h(k)^2) - k_f/(2*h(k)) ;
count = count + 1 ;
%j (diagonal) 
ivec(count) = j ; jvec(count) = j ; avec(count) = -2*r(j)*k_f/(h(k)^2) ;
count = count + 1 ;
%j+1 (right of the diagonal)
ivec(count) = j ; jvec(count) = j+1 ; avec(count) = r(j)*k_f/(h(k)^2) + k_f/(2*h(k)) ;
count = count + 1 ;
end

% Diriclet Boundary
ivec(count) = N ; jvec(count) = N; avec(count) = 1;

A = sparse(ivec,jvec,avec);
tol = 1e-4;
[T, counts] = GaussSeidel(A,bvec,tol);
GS_iters(k)=counts
% comparing analytical to numerical
errors(k) = norm(T - T_f)
figure(k); plot(r,T,'bo',r,T_f,'r')
xlabel('Radius in Pellet [m]')
ylabel('Temperature K')
title(['Comparing Numerical and Analytical Solution for ',num2str(N),' grid points'])
legend('numerical','analytical')
end
% number of iterations vs N grid points
figure(k+1);plot(Ns,GS_iters,'r-o')
xlabel('Number of Grid Points')
ylabel('Number of Iterations in GS Algorithm)')
tilte('Comparing GS Iterations Required vs Number of Grid Points')

% showing quadratic convergence for central FD between grid spacing and
% error
figure(k+2);plot(log(h),log(errors),'k-o')
xlabel('log(h)')
ylabel('log(errors(h))')
title('Log-Log Plot for Grid Spacing vs Error')

% Linear Polyfit to try and see slope of previous plot, slope =2 means
% quadratic convergence is true
coefs = polyfit(log(h),log(errors),1)
