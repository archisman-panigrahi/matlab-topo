% A point is represented by coordinates (i,k,l)
% In 2D, (i,k) is mapped to (k-1)*m + i, where k varies from 1 to n, and i from 1 to m
% In 3D, (i,k,l) is mapped to (l-1)*m*n + (k-1) * m + i, where l varies from 1 to p, k from 1 to n, i from 1 to m
% L_x = m, L_y = n, L_z = p 
%The symbol j is reserved for sqrt(-1)

% check whether lattice size is preassined. Otherwise assign initial values as below.
if ~exist('m','var')
	m = 8;
end

%Show the values
m

%Initialize
CX1Dp = zeros(m,m); %1D periodic Cosine matrix 
CX1Dnp = zeros(m,m);
SX1Dp = zeros(m,m);
SX1Dnp = zeros(m,m);

%1D matrices
for i = 1:m-1
	CX1Dp(i+1 , i) = 1/2;
	CX1Dp(i, i+1) = 1/2;

	SX1Dp(i+1, i) = j/2;
	SX1Dp(i, i+1) = -j/2;

	CX1Dnp(i+1 , i) = 1/2;
	CX1Dnp(i, i+1) = 1/2;

	SX1Dnp(i+1, i) = j/2;
	SX1Dnp(i, i+1) = -j/2;
endfor

% For periodic systems
	CX1Dp(1 , m) = 1/2;
	CX1Dp(m, 1) = 1/2;

	SX1Dp(1, m) = j/2;
	SX1Dp(m, 1) = -j/2;

% Constant matrix
M1D = eye(m);

