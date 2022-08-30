function fi=limiter_vanAlbada(q,r)

fi=1.;
if q > 0.
    fi=1.-q*(1.+q)*(1.-r)/(1.+q*q);
end
return;