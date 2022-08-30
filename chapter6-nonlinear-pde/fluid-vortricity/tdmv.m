function v=tdmv(a,b,c,x)

nx=length(x);
v(1)=0.;
for n=2:nx-1
    v(n)=a(n)*x(n-1)+b(n)*x(n)+c(n)*x(n+1);
end
v(nx)=0.;
return;