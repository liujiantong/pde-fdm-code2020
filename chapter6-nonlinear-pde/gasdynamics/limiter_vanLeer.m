function fi=limiter_vanLeer(q,r)

fi=1.;
if q > 0.
    fi=1.-2.*q*(1.-r)/(1.+q);
end
return;