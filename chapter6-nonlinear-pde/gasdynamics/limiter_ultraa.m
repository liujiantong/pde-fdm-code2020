function fi=limiter_ultraa(q,r)

fi=1.;
if q > 0. & q < r/(1.-r)
    fi=1.-2.*q*(1.-r)/r; return;
end
if q >= r/(1.-r)
    fi=-1.;
end
return;