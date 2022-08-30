function x = trid(A,b)
% solve tridiagonal system of equations
N = size(A,2);
for m = 2:N % Upper Triangularization
tmp = A(m,m - 1)/A(m - 1,m - 1);
A(m,m) = A(m,m) -A(m - 1,m)*tmp; A(m,m - 1) = 0;
b(m,:) = b(m,:) -b(m - 1,:)*tmp;
end
x(N,:) = b(N,:)/A(N,N);
for m = N - 1: -1: 1 % Back Substitution
x(m,:) = (b(m,:) -A(m,m + 1)*x(m + 1))/A(m,m);
end

