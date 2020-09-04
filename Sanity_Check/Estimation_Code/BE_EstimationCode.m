%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           BEHAVIOURAL EPIDEMIOLOGY
%                               ESTIMATION CODE
% Last Edited: 21st August, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETERS OF THE MODEL - ESTIMATED FROM DATA + CALIBRATED (reference: fwest.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COST PARAMETERS
%%%%%%%%%%%%%%%%
%pp = Phi_+        -> instantenous costs of infecting another person
%pm = Phi_-        -> instantenous costs of getting infected
%c =               -> flow costs of social distancing
%kap = Kappa       -> tuning parameter to smooth the agent's social distancing problem

%AGENT CONTACT RATES
%%%%%%%%%%%%%%%%%%%%
%bw = Beta_w       -> contact rate in productive activities
%bs = Beta_s       -> contact rate in social activities

%INFECTION/HOSPITALIZATION/TEST PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ti =              -> average length of a mild infection
%th =              -> average length of hospital stay
%mi =              -> hospitalization rate
%mh =              -> death rate. ifr rate initially set to 0.34%
%X = test capacity -> X(t) is the no. of tests available at time t. Maximum capacity: 1.000,000 per day.
%g = Gamma         -> efficiency tracing

%VACCINATION PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%
%p =               -> Pdf of vaccine arrival time. p(t) = probability vaccine arrives at t
%q =               -> Prob the vaccine has not arrived until time t. q(t) =sum_{k >= t}p(k)

%FIXED VALUES
%%%%%%%%%%%%%
%e0                -> fraction of people infected at the outset
%Lmin              -> Lower bound on Lockdown (essential services must remain open).
%da = Delta_A      -> agents' discount factor
%dp = Delta_P      -> government's discount factor
%T =               -> No. of time-periods in the model

%AUXILLIARY PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%
%Tm                -> expected time when the vaccine will arrive
%Tv                -> variance of the vaccine arrival time
%pop               -> size of population
%Ns                -> time before which the vaccine will not arrive for sure, i.e. the vaccine arrival time is supported on t:t=Ns,...,T
%Nt                -> time after which the vaccine will arrive with probability less than 1e-4 OR vaccine will arrive with prob almost 1 before Nt

%DATA VARIABLES
%%%%%%%%%%%%%%%
%POSCUMr           -> Cumulative no. of positive tests at time t
%NEGCUMr           -> Cumulative no. of negative tests at time t
%DCUMr             -> Cumulative no. of deaths at time t
%Hr                -> Daily no. of hospitalizations at time t
%POSNr             -> Daily no. of positive tests at time t (RAW)
%NEGNr             -> Daily no. of negative tests at time t (RAW)
%DNr               -> Daily no. of death at time t (RAW)
%POSN              -> Daily no. of positive tests at time t (SMOOTHED)
%NEGN              -> Daily no. of negatice tests at time t (SMOOTHED)
%DN                -> Daily no. of deaths at time t (SMOOTHED)
%HST               -> Daily no. of hospitalizations are time t (SMOOTHED)
%Tr                -> time before there was almost no testing, i.e. the number of negative tests was exactly zero
%Tmax              -> total number of data points in the dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CODE BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
global ti th mi mh e0 da kap Lmin rl1 beta_share 
global T Tmax p q X pop POSN NEGN DN HST Weight

%PATHNAME REQUIRED - change as necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_EstimationPlots = '../Figures/Simulation_Plots/';
path_SetupParameters  = '../Parameters/BE_DataSetup/';
path_EstimationParameters = '../Parameters/BE_EstimationCode/';
path_OldParameters = '../Parameters/Old_Parameters/';
path_SetupCode = '../Setup_Code/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(fullfile(path_SetupParameters, 'Initial_Parameters.mat'), 'file')==0
    run(fullfile(path_SetupCode, 'BE_DataSetup.m'));
    load(fullfile(path_SetupParameters, 'Initial_Parameters.mat'));
else
    load(fullfile(path_SetupParameters, 'Initial_Parameters.mat'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VARIABLE ORDER - Beta, Gamma, Phi_+, Phi_-/Phi_+, rl2, c
Pars_LB  = [0.200, 0.050, 0.000, 0.010, 0.990, 0.010];
Pars_UB = [0.400, 0.250, 200.000, 0.010, 0.99, 0.050];
% Pars_LB = [0.5  0.20    0.05   0.0     0.1     5    0.95    100.   0.01   ti th mi mh];
% Pars_UB = [0.5  0.40    0.25   200.0   0.1    5   0.995    100.   0.90   ti th mi mh];
NoOfIterations = 1;
Optimal_Parameters = zeros(length(Pars_UB), NoOfIterations);
Value_Function = zeros(NoOfIterations, 1);

for i=1:NoOfIterations
    Options = optimoptions('particleswarm','SwarmSize',100,'InitialSwarm',[],'UseVectorized',true,'Display','Off','MaxIterations',50,'MaxStallIterations',20,'FunctionTolerance',1e-12,'PlotFcn',@pswplotbestf);   
    Optimal_Parameters(:,i) = particleswarm(@(pars) Objective_Estimation(pars), length(Pars_LB), Pars_LB, Pars_UB, Options);  
    disp(Optimal_Parameters(:,i)) 
    [SS,XP0,XN0,D,SA,SR,L,tau,alpha] = Objective_Estimation(Optimal_Parameters(:,i)');
    Value_Function(:,i)= SS;
    figure
    subplot(2,2,1)
    plot((XP0(1:Tmax))*pop*1e6,'r','LineWidth',2.5);
    hold on
    plot((POSN(1:Tmax)),'b','LineWidth',2.5)
    hold off
    subplot(2,2,2)
    plot((D(1:Tmax))*pop*1e6,'r','LineWidth',2.5);
    hold on
    plot((DN(1:Tmax)),'b','LineWidth',2.5)
    hold off
    subplot(2,2,3)
    plot(SR(1:Tmax)./SA(1:Tmax),'r','LineWidth',2.5);
    hold on
    plot(L(1:Tmax),'b','LineWidth',2.5)
    hold off
    subplot(2,2,4)
    A = reshape(alpha,T,T);
    mesh(A(1:Tmax,1:Tmax))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









