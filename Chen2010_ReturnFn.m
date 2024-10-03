function F=Chen2010_ReturnFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)

F=-Inf;

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
    % Chen has typo in eqns, missing Tr (only gives transfers to homeowners), both here and in his codes Tr goes to both renters and homeowners.
    % Budget constraint
    cspend=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime; % -hprime=0, so Chen (2010) omits it, but I leave it here
    % cspend=c+p*d (consumption goods plus housing services)
    % Analytically, we can derive the split of cspend into c and p*d as
    c=cspend/(1+(p^(upsilon/(upsilon-1)))*((theta/(1-theta))^(1/(upsilon-1))));
    d=(cspend-c)/p;
    if c>0
        % utility function
        if upsilon==0
            uinner=(c^theta)*(d^(1-theta)); 
            % Chen2010 actually uses this, the formula in paper is kind of misleading as it just errors for upsilon=0 which is what is used. 
            % CES utility with unit elasticity of substitution (with upsilon=0) simplies to Cobb-Douglas utility.
        else
            uinner=(theta*(c^upsilon)+(1-theta)*(d^upsilon))^(1/upsilon);
        end
        F=(uinner^(1-sigma))/(1-sigma);
    end
    % Impose borrowing constraints
    if aprime<0
        F=-Inf;
    end
%% Owner
elseif hprime>0
    % Budget constraint
    c=(1+r)*a+(1-tau_p)*earnings+(1-delta_o)*h-tau_hhprime+(agej>=Jr)*b+Tr-aprime-hprime;
    if c>0
        % utility function
        if upsilon==0
            uinner=(c^theta)*(hprime^(1-theta)); 
            % Chen2010 actually uses this, the formula in paper is kind of misleading as it just errors for upsilon=0 which is what is used. 
            % CES utility with unit elasticity of substitution (with upsilon=0) simplies to Cobb-Douglas utility.
        else
            uinner=(theta*(c^upsilon)+(1-theta)*(hprime^upsilon))^(1/upsilon);
        end
        F=(uinner^(1-sigma))/(1-sigma);
    end
    % Impose collateral constraints
    if aprime<-(1-gamma)*hprime % 'mortgage' cannot be greater than (1-gamma)*hprime (can be thought of as the downpayment requirement)
        F=-Inf;
    end
end
