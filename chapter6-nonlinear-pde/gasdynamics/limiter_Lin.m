function fi=limiter_Lin(q,r)

fi=1.;
if q > 0. & q < 0.5
    fi=1.-2.*q*(1.-r); return;
end
if q >= 0.5 & q <= 1.
    fi=r; return;
end
if q > 1.
    fi=-1.+(1.+r)*exp(-50.*(1.-r)*(q-1.)); return;
end
return;