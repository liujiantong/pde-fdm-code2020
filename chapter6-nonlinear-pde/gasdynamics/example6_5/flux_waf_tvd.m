% Solves Riemann problem and computes fluxes for method WAF (TVD)

% r1, u1, p1   parameters in left cell
% r2, u2, p2   parameters in right cell
% rf, uf, pf   parameters on the interface
% c1, c2       sound speeds in left and right cells
% s1, sr       velocties of the fronts of S(hock)W(ave)s
% shl, shr     velocties of the fronts of the heads of R(arefaction)Ws
% stl, str     velocties of the tails of RWs

function [rw,uw,pw,cw]=flux_waf_tvd(sl1,sl2,sl3,sr1,sr2,sr3,g)

r1=sl1; u1=sl2/r1; p1=(1./g(8)-1.)*(sl3-0.5*sl2*u1);
r2=sr1; u2=sr2/r2; p2=(1./g(8)-1.)*(sr3-0.5*sr2*u2);

% sound speeds in left and right cells
c1 = sqrt(p1/(r1*g(8)));
c2 = sqrt(p2/(r2*g(8)));
du=u2-u1;

% Vacuum solution
uvac=g(4)*(c1+c2)-du;
if uvac <= 0.
    disp(' Vacuum generated by given data');
    return;
end    

% Apply Newton method
% Initial guesses
pv=0.5*(p1+p2)-0.125*du*(r1+r2)*(c1+c2);
pmin=min(p1,p2); pmax=max(p1,p2); qrat=pmax/pmin;
if ( qrat <= 2. ) & ( pmin <= pv & pv <= pmax )
    p=max(1.0e-6,pv);
else
    if pv < pmin
        pnu=c1+c2-g(7)*du;
        pde=c1/(p1^g(1))+c2/(p2^g(1));
        p=(pnu/pde)^g(3);
    else
        gel=sqrt((g(5)/r1)/(g(6)*p1+max(1.0e-6,pv)));
        ger=sqrt((g(5)/r2)/(g(6)*p2+max(1.0e-6,pv)));
        p=(gel*p1+ger*p2-du)/(gel+ger); p=max(1.0e-6,p);
    end
end    
p0=p; err=1.;
while err > 1.0e-5
    if p <= p1
        prat=p/p1; fl=g(4)*c1*(prat^g(1)-1.); ffl=1./(r1*c1*prat^g(2));
    else
        a=g(5)/r1; b=g(6)*p1; qrt=sqrt(a/(b+p));
        fl=(p-p1)*qrt; ffl=(1.-0.5*(p-p1)/(b+p))*qrt;
    end
    if p <= p2
        prat=p/p2; fr=g(4)*c2*(prat^g(1)-1.); ffr=1./(r2*c2*prat^g(2));
    else
        a=g(5)/r2; b=g(6)*p2; qrt=sqrt(a/(b+p));
        fr=(p-p2)*qrt; ffr=(1.-0.5*(p-p2)/(b+p))*qrt;
    end
    p=p-(fl+fr+du)/(ffl+ffr);
    err=2.*abs(p-p0)/(p+p0); p0=p;
end
if p0 < 0.
    p0=1.0e-6;
end

% Velocity of CD 
u0=0.5*(u1+u2+fr-fl);

uw(1)=u1; uw(4)=u2;
pw(1)=p1; pw(4)=p2;
rw(1)=r1; rw(4)=r2;
cw(2)=u0;
if p0 < p1
    shl=u1-c1; cw(1)=shl;
    cml=c1*(p0/p1)^g(1); stl=u0-cml;
    if stl <= 0.
        rw(2)=r1*(p0/p1)^g(8); uw(2)=u0; pw(2)=p0;
    else
        c=g(5)*(c1+g(7)*u1);
        uw(2)=c;
        rw(2)=r1*(c/c1)^g(4);
        pw(2)=p1*(c/c1)^g(3);
    end
else
    pml=p0/p1; sl=u1-c1*sqrt(g(2)*pml+g(1)); cw(1)=sl;
    rw(2)=r1*(pml+g(6))/(pml*g(6)+1.); uw(2)=u0; pw(2)=p0;
end
if p0 >= p2
    pmr=p0/p2; sr=u2+c2*sqrt(g(2)*pmr+g(1)); cw(3)=sr;
    rw(3)=r2*(pmr+g(6))/(pmr*g(6)+1.); uw(3)=u0; pw(3)=p0;
else
    shr=u2+c2; cw(3)=shr;
    cmr=c2*(p0/p2)^g(1); str=u0+cmr;
    if str >= 0.
        rw(3)=r2*(p0/p2)^g(8); uw(3)=u0; pw(3)=p0;
    else
        c=g(5)*(c2-g(7)*u2);
        uw(3)=g(5)*(-c2+g(7)*u2); rw(3)=r2*(c/c2)^g(4); pw(3)=p2*(c/c2)^g(3);
    end
end
return;