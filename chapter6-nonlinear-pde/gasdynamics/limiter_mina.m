function fi=limiter_mina(q,r)

fi=1.;
if q > 0. & q < 1.
    fi=1.-q*(1.-r); return;
end
if q >= 1.
    fi=r; return;
end
return;