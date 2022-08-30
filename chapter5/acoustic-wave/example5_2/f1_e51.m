function y=f1_e51(t)

p=2;
f0=(1+p)*((1+p)/p)^p;
y=0.;
if t < 1.
   y=f0*t*(1.-t)^p;
end
return;