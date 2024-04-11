m = 2, n = 3

xmatrix2D = diag(kron(kron(linspace(1,1,n), linspace(1,m,m)), [1, 1]))
ymatrix2D = diag(kron(kron(linspace(1,n,n), linspace(1,1,m)), [1, 1]))
