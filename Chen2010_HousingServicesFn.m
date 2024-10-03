function hs=Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)

hs=0;

% price of housing
p=(r+delta_r)/(1+r);
% wage
w=(1-alpha)*((r+delta_k)/alpha)^(alpha/(alpha-1));

% housing transactions costs
if hprime==h
    tau_hhprime=0;
else
    tau_hhprime=phi*h;
end

% earnings
earnings=w*kappaj*z;

%% Renter
if hprime==0
    % Budget constraint
    cspend=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime; % -hprime=0, so Chen (2010) omits it, but I leave it here [Chen has typo in eqns, missing Tr, it is here and in his codes]
    % cspend=c+p*d (consumption goods plus housing services)
    % Analytically, we can derive the split of cspend into c and p*d as
    c=cspend/(1+(p^(upsilon/(upsilon-1)))*((theta/(1-theta))^(1/(upsilon-1))));
    d=(cspend-c)/p;

    hs=d;
%% Owner
elseif hprime>0
    % Budget constraint
    % c=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime;

    hs=hprime;
end
