% Chen (2010) - A life-cycle analysis of social security with housing
%
% Two endogenous states: a (financial assets), h (housing)
% One exogenous state: z (stochastic earnings)
%
% Note that the decision variable d, the amount of housing-services
% purchased by a renter, can be solved for analytically and so does not
% need to be included in the code as a decision variable (instead the
% analytic solution is just included in the return function as a formula)

% Model has deterministic productivity growth g=0.015. You should
% renormalize the model (to solve of the balanced growth path), but I can't
% be bothered so I just ignore it (essentially I solve as if g=0).

% I renamed two things: 
%   Chen calls the stochastic process on earnings eta, I call it z.
%   Chen calls the deterministic earnings as function of age epsilonj, I call it kappaj

% The description of pensions in Chen (2010) is a bit of a mess.
% The budget constraints include a term I(j)*b, which is presumably the pension.
% There is a tau_p, which is presumably the payroll-tax.
% Page 604 says: "In the initial steady state, we choose the replacement rate vartheta so that the payroll-tax 
%         rate matches it empirical counterpart. Currently, the OASI (Old-Age and Survivors Insurance) 
%         rate is 10.7 percent.22 This implies vartheta = 0.483."
% I interpret this to mean that tau_p=0.107, and set b to balance the pension budget in eqm.
% Chen (2010) relates b to vartheta in Definition 1 in appendix, confirming the above.
% I eliminate vartheta from the model as it is not needed for anything [Having seen codes of Chen2010, he sometimes fixed
% tau_p and solved for b, but sometimes fixed vartheta and solved for tau_p (and b) in eqm. This didn't make final paper, hence 
% why he has vartheta which is essentially redundant.].

% Chen interpolates choice of assets, but forces choice of housing onto
% grid (his Appendix A.2). Here both are on grid, but grid for assets is vastly larger than
% that used by Chen. We use exact same grid on housing as Chen codes, 35 points equally spaced from 0 to 20.


%% 
n_d=0; % no decision variables
n_a=[301,35]; % Assets, Housing % Note: n_a(2)=35 is what Chen used
n_z=7; % Chen (2010) uses 7 points and Tauchen. I use Kirkby-Farmer-Tanaka-Toda (for discretizing life-cycle AR(1)) but stick to 7 points (even though this is less than I would like).

N_j=66; % total number of periods
Params.agejshifter=19; % starting age 20
% So model is annual for ages 20 to 85
Params.J=N_j;

useTauchen=1; % 0 or 1
% Chen2010 did discretization with Tauchen method, I prefer to use KFTT which is better

%% Parameters

% Discount factor
Params.beta=0.9622;

% Preferences
Params.sigma=2; % curvature of utility
Params.upsilon=0; % =0 implies a unit elasticity of substitution between the two types of consumption (consumption good and housing services)
Params.theta=0.8954; % share of nondurable consumption in the utility function

% Demographics
Params.Jr=46; % retire at age 65
Params.agej=1:1:Params.J; % agej is the period
Params.n=0.01; % Population growth rate

% Deterministic productivity growth
Params.g=0.015;

% Production function
Params.alpha=0.2743; % share of capital in Cobb-Douglas production fn
% Depreciation rates
Params.delta_r=0.0254; % annual depreciation rate for rental housing
Params.delta_o=0.013; % annual depreciation rate for owner-occupied housing
Params.delta_k=0.0951; % annual depreciation rate for physical capital (financial assets)

% Housing
Params.phi=0.05; % transaction cost for changing amount of housing
Params.gamma=0.2; % downpayment ratio (when buying a house)

% Earnings
Params.rho_z=0.96; % autocorrelation coeff of AR(1) on earnings
Params.sigma_z_epsilon=sqrt(0.045); % standard deviation of innovations to the AR(1) on earnings
Params.sigma_z1=sqrt(0.38); % standard deviation of initial earnings (at age 20, period 1)
% Note: asymptotic std dev of z is sqrt((0.045^2)/(1-0.96^2))=0.1607
% Chen (2010) codes make it clear that 0.045 is the variance of the innovations to z [paper says it is the std dev]
% Chen (2010) codes make it clear that 0.38 is the variance of normal initial dist on log(z) [paper omits the log()]

% Pension
Params.vartheta=0.483; % Not used for anything (Chen2010 reports this a the 'replacement rate', but is not used here)

% Taxes
Params.tau_p=0.107;

%First the raw data of Hansen 1993:
%Table II. Weights assigned to age-sex
%         Males       Females
%Age      Weight      Weight
%16-19    0.56        0.52
%20-24    0.78        0.69
%25-34    1.14        0.89
%35-44    1.37        0.90
%45-54    1.39        0.87
%55-64    1.33        0.84
%65 +     0.89        0.66
%14-17    0.56        0.52
%14-19    0.56        0.52
%18-24    0.78        0.69
%25-44    1.24        0.89
%45-64    1.37        0.86
% kappaj is set based on the first entries of the data for males
% Rather than use this, Chen (2010) codes "./benchmark/demographics.f90"
% contains "Age-efficiency Units from Hansen (1991)" which I copy here
Params.kappaj=[1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527, 1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, 1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, 1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658, 1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, 1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, 1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, 1.7701, 1.7048, 1.6396];
Params.kappaj=[Params.kappaj,zeros(1,Params.J-Params.Jr+1)]; % zeros for retirement

% Survival probabilities
% F. C. Bell and M. L. Miller (2005), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 120, Office of the Chief Actuary
% http://www.socialsecurity.gov/oact/NOTES/s2000s.html
% Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
% The raw data from there is
%          Sex and Exact Age
%     |  Male                                                Female
% Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
% 2010| [0.00587,0.00116,0.01086,0.01753,0.02785,0.39134]    [0.00495,0.00060,0.00734,0.01201,0.01912,0.34031]
% I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');
Params.sj=1-dj_temp(20:85);
Params.sj(1)=1;
Params.sj(end)=0;
% I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% I have additionally imposed that the prob of death at age 20 be zero and that prob of death at age 85 is one.
% Chen (2010) replication codes contains a "./benchmark/surv.txt" contains
% the survival probabilities he used, they are quite different, so I now provide those
Params.sj=[0.9985995, 0.9985575, 0.998498, 0.9985075, 0.9986075, 0.998641, 0.998743, 0.9987235, 0.9987385, 0.998705, 0.9986825, 0.9986935, 0.9986795, 0.9986615, 0.998573, 0.9984555, 0.998305, 0.9981125, 0.9979675, 0.997858, 0.9977755, 0.997652, 0.997384, 0.99714, 0.9969315, 0.9967195, 0.99666, 0.996468, 0.9962415, 0.9959035, 0.995637, 0.9953435, 0.9950225, 0.9946615, 0.9942525, 0.9937795, 0.993245, 0.9926615, 0.992031, 0.9913415, 0.9906035, 0.9897745, 0.9887855, 0.987602, 0.9862635, 0.9847635, 0.9832185, 0.981763, 0.9804655, 0.979231, 0.977836, 0.976205, 0.974443, 0.9725355, 0.9704165, 0.967929, 0.9650565, 0.9619165, 0.9585285, 0.954784, 0.9506335, 0.945818, 0.940017, 0.933017, 0.9249375, 0.9160355];
Params.sj(end)=0; % I need this for the accidental bequests calculations. Not sure how Chen (2010) treated it.

%% Grids and exogenous shocks
% Grid for housing, from 0 to 20 (which are min and max values used by Chen2010 codes)
minh=0; % min for housing; note that h'=0 means renter
maxh=20; % maximum value of housing: 
h_grid=minh + (maxh-minh)*linspace(0,1,n_a(2))'.^3; % note: ^3 is adding curvature, so most points are near hmin

% Note: -(1-Params.gamma)*maxh is the minimum possible value of assets, and maxassets is the maximum possible value
minassets=-(1-Params.gamma)*maxh;
maxassets=100; % maximum value of assets (Chen 2010 codes uses 100)
% First, the negative part of the grid: equally spaced from -(1-Params.gamma)*maxh up to 0. 
n_a_neg = round(0.1*n_a(1)); % use 1/10 of points in the negative assets range
asset_grid_neg = linspace(minassets,0,n_a_neg)';
% Second, the positive part of the grid (more points nearer 0)
asset_grid_pos=maxassets*linspace(0,1,n_a(1)-n_a_neg+1)'.^3;
asset_grid = [asset_grid_neg(1:n_a_neg-1); asset_grid_pos]; % Note: both contain zero, so omit it from asset_grid_neg
a_grid=[asset_grid; h_grid]; % stacked column vector

if useTauchen==0
    % Chen2010 used Tauchen method, we instead use KFTT which will be better
    kfttoptions.initialj1sigmaz=Params.sigma_z1;
    kfttoptions.nSigmas=sqrt(7/2); % plus/minus 7/2 std deviations as the max/min grid points % This is what Chen2010 codes use (in Tauchen)
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs] = discretizeLifeCycleAR1_KFTT(0,Params.rho_z,Params.sigma_z_epsilon,n_z,N_j,kfttoptions);
    z_grid_J=exp(z_grid_J);
elseif useTauchen==1
    % Just for comparison, use Tauchen to see what the grid Chen (2010) would have had looks like
    tauchenoptions=struct();
    Tauchen_q=sqrt(7/2); % plus/minus sqrt(7/2) std deviations as the max/min grid points % This is exactly what Chen (2010) codes use
    [z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho_z,Params.sigma_z_epsilon,n_z,Tauchen_q,tauchenoptions);
    z_grid=exp(z_grid); % Ranges from 0.7 to 1.4
    [z_mean,~,~,~]=MarkovChainMoments(z_grid,pi_z);
    z_grid=z_grid/z_mean; % renormalize grid to be mean 1 (Chen2010 codes do this)
    % I will put the initial dist onto this grid
    jequaloneDistz=MVNormal_ProbabilitiesOnGrid(log(z_grid),1,Params.sigma_z1,n_z); % Chen2010 codes do initial dist on z as log-normal
end

d_grid=[]; % no d variable

%% Return fn
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)...
    Chen2010_ReturnFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);

%% General eqm parameters
Params.r=0.0635; % This is an initial guess, but it is also the value that Chen (2010) reports as being the initial stationary general eqm
Params.Tr=0.05; % accidental bequests of assets and housing

% Three other parameters are being determined in general equilibrium, but
% we can shortcut for these. They are b, p and w.

% p and w are both determined in general eqm. But since both are just
% simple functions of other general eqm parameters, we can just hard-code
% them (inside return function and elsewhere).
Params.p=(Params.r+Params.delta_r)/(1+Params.r); % price of housing
Params.w=(1-Params.alpha)*((Params.r+Params.delta_k)/Params.alpha)^(Params.alpha/(Params.alpha-1)); % wage

% Because labor supply is exogenous you can actually figure out of b is without having to
% solve general eqm (point 4 in Definition 1 of his Appendix A.1).
Params.b=0.1;
% This is an initial guess, and then we set b to clear the general eqm
% below. It is just below the first creation of AggVars.


%% Solve the value function
vfoptions.divideandconquer=1;
vfoptions.level1n=[9,n_a(2)];
tic;
if useTauchen==0
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
elseif useTauchen==1
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
end
vftime=toc


%% Age distribution
AgeWeightParamNames={'mewj'};
Params.mewj=ones(1,N_j);
for jj=1:N_j-1
    Params.mewj(jj+1)=Params.mewj(jj)*Params.sj(jj)/(1+Params.n); % iteratively apply formula that mewj is determined by conditional survival probabilities (sj) and population growth rate (n)
end
Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one

%% Initial age 20 (period 1) distribution
jequaloneDist=zeros([n_a,n_z],'gpuArray');
% Everyone is born with zero assets and zero housing
[~,zeroassetindex]=min(abs(asset_grid));
jequaloneDist(zeroassetindex,1,:)=shiftdim(jequaloneDistz,-2); % initial dist of z, with zero assets and zero housing

%% Agent distribution
simoptions=struct(); % defaults
if useTauchen==0
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
elseif useTauchen==1
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
end

%% AggVars
FnsToEvaluate.A=@(aprime,hprime,a,h,z) a;
FnsToEvaluate.N=@(aprime,hprime,a,h,z,kappaj) kappaj*z;
FnsToEvaluate.H=@(aprime,hprime,a,h,z) h; % housing
FnsToEvaluate.Hr=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr) (hprime==0)*Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr); % rental housing=housing services used by renters
FnsToEvaluate.pensiontaxrevenue=@(aprime,hprime,a,h,z,tau_p,w,kappaj) tau_p*w*kappaj*z;
FnsToEvaluate.pensionspend=@(aprime,hprime,a,h,z,agej,Jr,b) (agej>=Jr)*b;
FnsToEvaluate.accidentalbeqleft=@(aprime,hprime,a,h,z,r,sj,delta_o) (1+r)*aprime*(1-sj)+(1-delta_o)*hprime*(1-sj);

if useTauchen==0
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
elseif useTauchen==1
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end

% Note: Because labor supply is exogenous we can set b to clear the general eqm constraint directly, rather than solving for it as part of general eqm
Params.b=Params.b*(AggVars.pensiontaxrevenue.Mean/AggVars.pensionspend.Mean);
% This is essentially what Chen (2010) does but he does it on the replacement rate (vartheta) rather than the pension amount (b).

%% Everything is working, now for general equilibrium
GEPriceParamNames={'r','Tr'};

GeneralEqmEqns.capitalmarkets=@(r,A,N,alpha,delta_k,Hr,delta_r) r-(alpha*((A-Hr*(1-((r+delta_r)/(1+r))))^(alpha-1))*(N^(1-alpha))-delta_k); % r=marginal product of capital, minus depreciation; with K'=A'-Hr'*(1-p), and p=(r+delta_r)/(1+r);
GeneralEqmEqns.accidentalbequests=@(Tr,accidentalbeqleft,n) Tr-accidentalbeqleft/(1+n); % Eqn A.2 from Appendix of Chen2010

%% Alright, solve for the general eqm
heteroagentoptions.verbose=1; % defaults
if useTauchen==0
    [p_eqm,~,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z_J, d_grid, a_grid, z_grid_J, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
elseif useTauchen==1
    [p_eqm,~,GECondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end

% Update Params based on the general eqm
Params.r=p_eqm.r;
Params.Tr=p_eqm.Tr;

% What Chen (2010) got:
% Params.r=0.064; % This is the value that Chen (2010) reports as being the initial stationary general eqm
% Params.Tr=0.050; % accidental bequests of assets and housing
% Obviously, since the code here just pretends g=0, we won't get the same thing.

Params.p=(Params.r+Params.delta_r)/(1+Params.r); % price of housing
Params.w=(1-Params.alpha)*((Params.r+Params.delta_k)/Params.alpha)^(Params.alpha/(Params.alpha-1)); % wage

%% A little bit of output about the general eqm

% We will compute some stats conditional on being a homeowner
simoptions.conditionalrestrictions.Homeowners=@(aprime,hprime,a,h,z) (hprime>0);

% Some additional outputs
FnsToEvaluate.Homeownership=@(aprime,hprime,a,h,z) (hprime>0);
FnsToEvaluate.TotalWealth=@(aprime,hprime,a,h,z) a+h; % NOT SURE IF THIS IS CORRECT DEFINITION OF TOTAL WEALTH
FnsToEvaluate.LoanToValueRatio=@(aprime,hprime,a,h,z) (aprime<0)*abs(aprime/hprime);
FnsToEvaluate.earnings=@(aprime,hprime,a,h,z,w,kappaj) w*kappaj*z;
FnsToEvaluate.consumption=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr) Chen2010_ConsumptionFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);
FnsToEvaluate.housingservices=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr) Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);

if useTauchen==0
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
elseif useTauchen==1
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end

K=AggVars.A.Mean-AggVars.Hr.Mean*(1-Params.p); % physical capital
Y=(K^Params.alpha)*(AggVars.N.Mean^(1-Params.alpha)); % output

% Chen (2010) Table 2 on page 605 reports the following
fprintf('Quantitative properties of the benchmark economy \n')
% Comments at end of each line are values in the paper for comparison
fprintf('Targeted Variables \n')
fprintf('Payroll tax rate is %1.3f \n',Params.tau_p) % 0.107
fprintf('r is %1.2f%% \n',100*Params.r) % 6.35%
fprintf('K/Y is %1.3f \n',K/Y) % 1.729
fprintf('H/Y is %1.3f \n',AggVars.H.Mean/Y) % 1.075
fprintf('Homeownership rate is %2.1f%% \n', 100*AggVars.Homeownership.Mean) % 65.0%
fprintf('Nontargeted Variables \n')
% fprintf('pH/(C+pH) is %2.1f \% \n',100*Params.p*AggVars.H.Mean/()) % 10.7% I SKIPPED AS CANT BE BOTHERED CREATING CONSUMPTION FN
fprintf('Ho/(A+Ho) is %2.1f%% \n',100*AggVars.H.Mean/(AggVars.A.Mean+AggVars.H.Mean)) % 32.2%  % I think this is the correct calculation, not sure
fprintf('Gini for total wealth is %1.2f \n',AllStats.TotalWealth.Gini) % 0.73
fprintf('Gini for financial wealth is %1.2f \n',AllStats.A.Gini) % 0.93
fprintf('Gini for housing is %1.2f \n',AllStats.H.Gini) % 0.52 % NOT SURE HOW RENTAL HOUSING SERVICES ARE TREATED HERE, GUESSING JUST AS ZEROS?
fprintf('Mean loan-to-value ratio (for borrowers) is %2.1f \n',100*AllStats.Homeowners.LoanToValueRatio.Mean) % 0.93


% Chen2010, bottom of pg 603 says that the Gini coeff for labor income is 0.40
simoptions.agegroupings=[1,Params.Jr]; % working age and retirees
if useTauchen==0
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
elseif useTauchen==1
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end
AgeConditionalStats.earnings.Gini(1) % is 0.33, decent way below 0.4
AllStats.earnings.Gini % 0.45, so still not 0.4

% Quick look at the replacement rate
Params.vartheta % 0.48
Params.b/AgeConditionalStats.earnings.Mean(1) % 0.46, so is correct

% Plot some life-cycle profiles to see more about what is going on
simoptions=struct();
if useTauchen==0
    AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
elseif useTauchen==1
    AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end

% Note: Chen (2010) graphs some life-cycle profiles based on 5-year age bins. 
simoptions.agegroupings=1:5:N_j; % 5-period bins [numbers are the start of each bin]
if useTauchen==0
    AgeConditionalStats5=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
elseif useTauchen==1
    AgeConditionalStats5=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
end


figure(1)
subplot(3,1,1); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.earnings.Mean)
title('Earnings (age conditional mean of)')
subplot(3,1,2); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.A.Mean)
title('Financial Wealth (age conditional mean of)')
subplot(3,1,3); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.H.Mean)
title('Housing (age conditional mean of)')

figure(2)
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.consumption.Mean)
hold on
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.housingservices.Mean)
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.earnings.Mean)
hold off
legend('consumption','housing services', 'earnings')
title('Consumption, Housing Services and Earnings')

% Plot home-ownership based on 5-year age bins
figure(3)
plot(20:5:85,AgeConditionalStats5.Homeownership.Mean)
title('Home-ownership rate [Fig 1 of Chen 2010]')


