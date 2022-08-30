function fi=limiter_supera(q,r)

fi=1.;
if q > 0. & q < 0.5
    fi=1.-2.*q*(1.-r); return;
end
if q >= 0.5 & q <= 1.
    fi=r; return;
end
if q > 1. & q < 2.
    fi=1.-q*(1.-r); return;
end
if q >= 2.
    fi=2.*r-1.;
end
return;