function [x, num_it ] = GaussSeidel(A,b,tol)
% This function accepts a matrix A and a vector b which are related by
% Ax = b. It uses the Gauss Seidel algorithm to iterate until a solution is
% achieved

% Check A is square
[M , N] = size(A) ;
if(M~=N)
    disp('Matrix is not square!')
    disp(['Matrix size is ' , num2str(M) ,' by ',num2str(N)])
    x = 0;
    return;
end
[P , Q] = size(b);
% Check sizes of A and b are well posed
if(Q~=1 || P~=M)
    disp('System is not correct size!')
    disp(['Vector should have 1 column and it has ', num2str(Q)])
    disp(['Matrix and vector should have the same number of rows! ', newline, ...
          'There are ', num2str(M), ' rows in the matrix and ',num2str(P), ...
          ' rows in vector'])
    x=0;
    return;
end

% Previous checks ensure b is a row vector and that it has the same number
% of rows of A. If these are true, then b and x have the same size.
% Initial guess vector with zeros is acceptable for Gauss Seidel.
xnew = zeros(size(b)) ; % starting solution vector
e = tol + 1;
count = 0 ;
while(e>tol)
    % update xold before computing new xnew
    xold = xnew;
    % do iteratiotn over each position in the vector
    for m=1:M
        % reinitialize sums and 
        sum = 0; % term based on last iteration
        % loop for contributtion from thtis iteration
        for n = 1:m-1
            sum = sum + A(m,n)*xnew(n) ; 
        end
        
        % loop for contribution from last iteration
        for n = m+1:M
            sum = sum + A(m,n)*xold(n) ;
        end
        xnew(m) = (b(m) - sum ) / A(m,m) ;
    end
    % check norm between new and old solution, end will interation if e<tol
    e = norm(xnew-xold) ; 
    count = count + 1;
end
 
x = xnew;
num_it = count;
end