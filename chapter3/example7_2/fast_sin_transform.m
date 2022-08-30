% Fast sine transform on input matrix A of dimension n-by-m

function y=fast_sin_transform(a)
[n,m]=size(a); 
b=[zeros(1,m);a;zeros(n+1,m)]; 
b=imag(fft(b)); 
% In Matlab vectors and matrices are indexed starting at 1, not 0
y=b(2:n+1,:);