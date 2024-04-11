% A point is represented by coordinates (i,k,l)
% In 2D, (i,k) is mapped to (k-1)*m + i, where k varies from 1 to n, and i from 1 to m
% In 3D, (i,k,l) is mapped to (l-1)*m*n + (k-1) * m + i, where l varies from 1 to p, k from 1 to n, i from 1 to m
% L_x = m, L_y = n, L_z = p 
%The symbol j is reserved for sqrt(-1)

% check whether lattice size is preassined. Otherwise assign initial values as below.
%Show the values
m,n,p

%Initialize
CX3Dnp = zeros(m * n * p,m * n * p);
SX3Dnp = zeros(m * n * p,m * n * p); %3D non-periodic Sine matrix describing hopping in x direction

CY3Dnp = zeros(m * n * p,m * n * p);
SY3Dnp = zeros(m * n * p,m * n * p);

CZ3Dnp = zeros(m * n * p,m * n * p);
SZ3Dnp = zeros(m * n * p,m * n * p);

%3D matrices

%3Dx
for l= 1:p
	for k= 1:n
		for i= 1:m-1
			CX3Dnp((l-1)*m*n + (k-1)*m + i + 1, (l-1)*m*n + (k-1)*m + i) = 1/2;
			CX3Dnp((l-1)*m*n + (k-1)*m + i, (l-1)*m*n + (k-1)*m + i + 1) = 1/2;

			SX3Dnp((l-1)*m*n + (k-1)*m + i + 1, (l-1)*m*n + (k-1)*m + i) = j/2;
			SX3Dnp((l-1)*m*n + (k-1)*m + i, (l-1)*m*n + (k-1)*m + i + 1) = -j/2;
		endfor
	endfor
endfor

%3Dy
for i= 1:m
	for l= 1:p
		for k= 1:n-1
			CY3Dnp((l-1)*m*n + k*m + i, (l-1)*m*n + (k-1)*m +i) = 1/2;
			CY3Dnp((l-1)*m*n + (k-1)*m + i, (l-1)*m*n + k*m + i) = 1/2;

			SY3Dnp((l-1)*m*n + k*m + i, (l-1)*m*n + (k-1)*m +i) = j/2;
			SY3Dnp((l-1)*m*n + (k-1)*m + i, (l-1)*m*n + k*m + i) = -j/2;
		endfor
	endfor
endfor

%3Dz
for i = 1:m
	for k = 1:n
		for l = 1:p-1
			CZ3Dnp(l*m*n + (k-1)*m + i, (l-1)*m*n + (k-1)*m +i) = 1/2;
			CZ3Dnp((l-1)*m*n + (k-1)*m + i, l*m*n + (k-1)*m +i) = 1/2;

			SZ3Dnp(l*m*n + (k-1)*m + i, (l-1)*m*n + (k-1)*m +i) = j/2;
			SZ3Dnp((l-1)*m*n + (k-1)*m + i, l*m*n + (k-1)*m +i) = -j/2;
		endfor
		%For periodic systems
	endfor
endfor

% Constant matrix
M3D = eye(m*n*p);