function [ A ] = Create_Poisson_problem_A( N )

  % Create the archtypical matrix A for an N x N Poisson problem (5-point
  % stencil.
  sz = N*N;
  A = zeros( sz, sz );

  % Set the diagonal
  for i=1:sz
      A(i,i) = 4;
  end
  
  % Set the entries of the first sub and super diagonals
  for i=1:sz-1
      if mod(i, N) ~= 0
          A(i,i+1) = -1;
      end
  end
  
  for i = 2:sz
      j = i - 1;
      if mod(j, N) ~= 0
          A(i, j) = -1;
      end
  end
      
  
  % Set the other off-diagonal entries
  
  for i = 1:sz-N
      A(i,i+N) = -1;
  end
  
  for i = N+1:sz
      A(i, i-N) = -1;
  end
  
end





