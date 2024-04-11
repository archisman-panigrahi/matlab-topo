% A point is represented by coordinates (i,k,l)
% In 2D, (i,k) is mapped to (k-1)*m + i, where k varies from 1 to n, and i from 1 to m
% In 3D, (i,k,l) is mapped to (l-1)*m*n + (k-1) * m + i, where l varies from 1 to p, k from 1 to n, i from 1 to m
% L_x = m, L_y = n, L_z = p 
%The symbol j is reserved for sqrt(-1)

% check whether lattice size is preassined. Otherwise assign initial values as below.
% if ~exist('m','var')
% 	m = 8;
% end
% if ~exist('n','var')
% 	n = 8;
% end

%Show the values
m,n

%Initialize
CX2Dp = zeros(m * n,m * n);
CX2Dnp = zeros(m * n,m * n);
SX2Dp = zeros(m * n,m * n); %2D non-periodic Sine matrix describing hopping in x direction
SX2Dnp = zeros(m * n,m * n);


CY2Dp = zeros(m * n,m * n);
CY2Dnp = zeros(m * n,m * n);
SY2Dp = zeros(m * n,m * n);
SY2Dnp = zeros(m * n,m * n);

% 2D matrices

%2Dx
for k = 1:n
	for i = 1:m-1
		CX2Dp((k-1)*m + i + 1, (k-1)*m + i) = 1/2;
		CX2Dp((k-1)*m + i, (k-1)*m + i + 1) = 1/2;

		SX2Dp((k-1)*m + i + 1, (k-1)*m + i) = j/2;
		SX2Dp((k-1)*m + i, (k-1)*m + i + 1) = -j/2;

		CX2Dnp((k-1)*m + i + 1, (k-1)*m + i) = 1/2;
		CX2Dnp((k-1)*m + i, (k-1)*m + i + 1) = 1/2;

		SX2Dnp((k-1)*m + i + 1, (k-1)*m + i) = j/2;
		SX2Dnp((k-1)*m + i, (k-1)*m + i + 1) = -j/2;
	endfor
% For periodic systems
	CX2Dp((k-1)*m + 1, (k-1)*m + m) = 1/2;
	CX2Dp((k-1)*m + m, (k-1)*m + 1) = 1/2;

	SX2Dp((k-1)*m + 1, (k-1)*m + m) = j/2;
	SX2Dp((k-1)*m + m, (k-1)*m + 1) = -j/2;
endfor

%2Dy
for i = 1:m
	for k = 1:n-1
		CY2Dp(k*m + i, (k-1)*m + i) = 1/2;
		CY2Dp((k-1)*m + i, k*m + i) = 1/2;

		SY2Dp(k*m + i, (k-1)*m + i) = j/2;
		SY2Dp((k-1)*m + i, k*m + i) = -j/2;

		CY2Dnp(k*m + i, (k-1)*m + i) = 1/2;
		CY2Dnp((k-1)*m + i, k*m + i) = 1/2;

		SY2Dnp(k*m + i, (k-1)*m + i) = j/2;
		SY2Dnp((k-1)*m + i, k*m + i) = -j/2;
	endfor
% For periodic systems
	CY2Dp(i, (n-1)*m + i) = 1/2;
	CY2Dp((n-1)*m + i, i) = 1/2;

	SY2Dp(i, (n-1)*m + i) = j/2;
	SY2Dp((n-1)*m + i, i) = -j/2;
endfor

% Constant matrix
M2D = eye(m*n);

%sin(kx)*sin(ky) matrix in real space

SXnpYnp2D = zeros(m * n, m * n);
SXpYnp2D = zeros(m * n, m * n);
SXnpYp2D = zeros(m * n, m * n);
SXpYp2D = zeros(m * n, m * n);

for a = 1:(m-1)
	for b = 2:n
		SXnpYnp2D(a + (b-1) * m, a + 1 + (b - 2) * m) = 1/4;
		SXnpYnp2D(a + 1 + (b - 2) * m, a + (b-1) * m) = 1/4;
	endfor
endfor

for a = 1:(m-1)
	for b = 1:(n-1)
		SXnpYnp2D(a + (b-1) * m, a + 1 + b * m) = -1/4;
		SXnpYnp2D(a + 1 + b * m, a + (b-1) * m) = -1/4;
	endfor
endfor

%Periodic only in X
for b = 2:(n-1)
	SXpYnp2D(b*m, (b-2) * m + 1) = 1/4;
	SXpYnp2D((b-2) * m + 1, b*m) = 1/4;

	SXpYnp2D(b*m, b*m + 1) = -1/4;
	SXpYnp2D(b*m + 1, b*m) = -1/4;
endfor

SXpYnp2D(n * m , (n-2) * m + 1) = 1/4;
SXpYnp2D((n-2) * m + 1, n * m) = 1/4;

SXpYnp2D(m, m + 1) = -1/4;
SXpYnp2D(m + 1, m) = -1/4;

%Periodic only in Y
for a = 1:(m-1)
	SXnpYp2D(a + (n - 1) * m, a + 1) = -1/4;
	SXnpYp2D(a + 1, a + (n - 1) * m) = -1/4;

	SXnpYp2D(a, a + 1 + (n - 1) * m) = 1/4;
	SXnpYp2D(a + 1 + (n - 1) * m, a) = 1/4;
endfor

%Both periodic
%Corners
SXpYp2D(m*n, 1) = -1/4;
SXpYp2D(1,m*n) = -1/4;

SXpYp2D(m, 1 + (n - 1) * m) = 1/4;
SXpYp2D(1 + (n - 1) * m, m) = 1/4;

%Finally add up
SXpYp2D += SXnpYp2D + SXpYnp2D + SXnpYnp2D;
SXnpYp2D += SXnpYnp2D;
SXpYnp2D += SXnpYnp2D;

%%%%%%%% sin(2kx) matrix %%%%%%

S2X2Dnp = zeros(m*n, m*n);
S2X2Dp = zeros(m*n, m*n);

for b = 1:n
	for a = 1:(m-2)
		S2X2Dnp(a + (b-1)*m, a + (b-1)*m + 2) = -j/2;
		S2X2Dnp(a + (b-1)*m + 2, a + (b-1)*m) = j/2;
	endfor
	%periodic
	S2X2Dp (b*m - 1, 1 + (b-1)*m) = -j/2;
	S2X2Dp (1 + (b-1)*m, b*m - 1) = j/2;

	S2X2Dp(b*m, 2 + (b-1)*m) = -j/2;
	S2X2Dp(2 + (b-1)*m, b*m) = j/2;
endfor

%Finally add up
S2X2Dp += S2X2Dnp;


%%%%%%%% sin(kx)*cos(ky) matrix %%%%%%%%%

SXnpCYnp2D = zeros(m * n, m * n);
SXpCYnp2D = zeros(m * n, m * n);
SXnpCYp2D = zeros(m * n, m * n);
SXpCYp2D = zeros(m * n, m * n);

for a = 1:(m-1)
	for b = 2:n
		SXnpCYnp2D(a + (b-1) * m, a + 1 + (b - 2) * m) = -j/4;
		SXnpCYnp2D(a + 1 + (b - 2) * m, a + (b-1) * m) = j/4;
	endfor
endfor

for a = 1:(m-1)
	for b = 1:(n-1)
		SXnpCYnp2D(a + (b-1) * m, a + 1 + b * m) = -j/4;
		SXnpCYnp2D(a + 1 + b * m, a + (b-1) * m) = j/4;
	endfor
endfor

%Periodic only in X
for b = 2:(n-1)
	SXpCYnp2D(b*m, (b-2) * m + 1) = -j/4;
	SXpCYnp2D((b-2) * m + 1, b*m) = j/4;

	SXpCYnp2D(b*m, b*m + 1) = -j/4;
	SXpCYnp2D(b*m + 1, b*m) = j/4;
endfor

SXpCYnp2D(n * m , (n-2) * m + 1) = -j/4;
SXpCYnp2D((n-2) * m + 1, n * m) = j/4;

SXpCYnp2D(m, m + 1) = -j/4;
SXpCYnp2D(m + 1, m) = j/4;

%Periodic only in Y
for a = 1:(m-1)
	SXnpCYp2D(a + (n - 1) * m, a + 1) = -j/4;
	SXnpCYp2D(a + 1, a + (n - 1) * m) = j/4;

	SXnpCYp2D(a, a + 1 + (n - 1) * m) = -j/4;
	SXnpCYp2D(a + 1 + (n - 1) * m, a) = j/4;
endfor

%Both periodic
%Corners
SXpCYp2D(m*n, 1) = -j/4;
SXpCYp2D(1,m*n) = j/4;

SXpCYp2D(m, 1 + (n - 1) * m) = -j/4;
SXpCYp2D(1 + (n - 1) * m, m) = j/4;

%Finally add up
SXpCYp2D += SXnpCYp2D + SXpCYnp2D + SXnpCYnp2D;
SXnpCYp2D += SXnpCYnp2D;
SXpCYnp2D += SXnpCYnp2D;

%%%%%%%% sin(ky)*cos(kx) matrix %%%%%%%%%

CXnpSYnp2D = zeros(m * n, m * n);
CXpSYnp2D = zeros(m * n, m * n);
CXnpSYp2D = zeros(m * n, m * n);
CXpSYp2D = zeros(m * n, m * n);

for a = 1:(m-1)
	for b = 2:n
		CXnpSYnp2D(a + (b-1) * m, a + 1 + (b - 2) * m) = j/4;
		CXnpSYnp2D(a + 1 + (b - 2) * m, a + (b-1) * m) = -j/4;
	endfor
endfor

for a = 1:(m-1)
	for b = 1:(n-1)
		CXnpSYnp2D(a + (b-1) * m, a + 1 + b * m) = -j/4;
		CXnpSYnp2D(a + 1 + b * m, a + (b-1) * m) = j/4;
	endfor
endfor

%Periodic only in X
for b = 2:(n-1)
	CXpSYnp2D(b*m, (b-2) * m + 1) = j/4;
	CXpSYnp2D((b-2) * m + 1, b*m) = -j/4;

	CXpSYnp2D(b*m, b*m + 1) = -j/4;
	CXpSYnp2D(b*m + 1, b*m) = j/4;
endfor

CXpSYnp2D(n * m , (n-2) * m + 1) = j/4;
CXpSYnp2D((n-2) * m + 1, n * m) = -j/4;

CXpSYnp2D(m, m + 1) = j/4;
CXpSYnp2D(m + 1, m) = -j/4;

%Periodic only in Y
for a = 1:(m-1)
	CXnpSYp2D(a + (n - 1) * m, a + 1) = -j/4;
	CXnpSYp2D(a + 1, a + (n - 1) * m) = j/4;

	CXnpSYp2D(a, a + 1 + (n - 1) * m) = j/4;
	CXnpSYp2D(a + 1 + (n - 1) * m, a) = -j/4;
endfor

%Both periodic
%Corners
CXpSYp2D(m*n, 1) = -j/4;
CXpSYp2D(1,m*n) = j/4;

CXpSYp2D(m, 1 + (n - 1) * m) = j/4;
CXpSYp2D(1 + (n - 1) * m, m) = -j/4;

%Finally add up
CXpSYp2D += CXnpSYp2D + CXpSYnp2D + CXnpSYnp2D;
CXnpSYp2D += CXnpSYnp2D;
CXpSYnp2D += CXnpSYnp2D;

%%%%%%%% cos(kx)*cos(ky) matrix %%%%%%%%%

CXnpCYnp2D = zeros(m * n, m * n);
CXpCYnp2D = zeros(m * n, m * n);
CXnpCYp2D = zeros(m * n, m * n);
CXpCYp2D = zeros(m * n, m * n);

for a = 1:(m-1)
	for b = 2:n
		CXnpCYnp2D(a + (b-1) * m, a + 1 + (b - 2) * m) = 1/4;
		CXnpCYnp2D(a + 1 + (b - 2) * m, a + (b-1) * m) = 1/4;
	endfor
endfor

for a = 1:(m-1)
	for b = 1:(n-1)
		CXnpCYnp2D(a + (b-1) * m, a + 1 + b * m) = 1/4;
		CXnpCYnp2D(a + 1 + b * m, a + (b-1) * m) = 1/4;
	endfor
endfor

%Periodic only in X
for b = 2:(n-1)
	CXpCYnp2D(b*m, (b-2) * m + 1) = 1/4;
	CXpCYnp2D((b-2) * m + 1, b*m) = 1/4;

	CXpCYnp2D(b*m, b*m + 1) = 1/4;
	CXpCYnp2D(b*m + 1, b*m) = 1/4;
endfor

CXpCYnp2D(n * m , (n-2) * m + 1) = 1/4;
CXpCYnp2D((n-2) * m + 1, n * m) = 1/4;

CXpCYnp2D(m, m + 1) = 1/4;
CXpCYnp2D(m + 1, m) = 1/4;

%Periodic only in Y
for a = 1:(m-1)
	CXnpCYp2D(a + (n - 1) * m, a + 1) = 1/4;
	CXnpCYp2D(a + 1, a + (n - 1) * m) = 1/4;

	CXnpCYp2D(a, a + 1 + (n - 1) * m) = 1/4;
	CXnpCYp2D(a + 1 + (n - 1) * m, a) = 1/4;
endfor

%Both periodic
%Corners
CXpCYp2D(m*n, 1) = 1/4;
CXpCYp2D(1,m*n) = 1/4;

CXpCYp2D(m, 1 + (n - 1) * m) = 1/4;
CXpCYp2D(1 + (n - 1) * m, m) = 1/4;

%Finally add up
CXpCYp2D += CXnpCYp2D + CXpCYnp2D + CXnpCYnp2D;
CXnpCYp2D += CXnpCYnp2D;
CXpCYnp2D += CXnpCYnp2D;

%%%%%%% sin(2ky) %%%%%%%
S2Y2Dnp = zeros(m*n, m*n);
S2Y2Dp = zeros(m*n, m*n);

for a = 1:m
	for b = 1:(n-2)
		S2Y2Dnp(a + (b - 1)*m, a + (b + 1)*m) = -j/2;
		S2Y2Dnp(a + (b + 1)*m, a + (b - 1)*m) = j/2;
	endfor
	S2Y2Dp(a, a + (n-2)*m) = j/2;
	S2Y2Dp(a + (n-2)*m, a) = -j/2;

	S2Y2Dp(a + m, a + (n-1)*m) = j/2;
	S2Y2Dp(a + (n-1)*m, a + m) = -j/2;
endfor
%Finally add up
S2Y2Dp += S2X2Dnp;
