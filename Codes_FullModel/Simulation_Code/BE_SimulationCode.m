%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           BEHAVIOURAL EPIDEMIOLOGY
%                               SIMULATION CODE
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
global ti th mi mh e0 dp da kap Lmin ifr kap
global bw bs g pp pm c
global T Tmax pop p q X  

%PATHNAME REQUIRED - change as necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_Plots = '../Figures/Simulation_Plots';
path_SetupCode  = '../Setup_Code/';
path_EstimationCode = '../Estimation_Code/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PARAMETER VALUES
%%%%%%%%%%%%%%%%%
%NOTE: Parameters set for ifr=0.0034, Lmin=0.3

% if exist(fullfile(path_EstimationCode, 'EstimationParameters.mat')) == 0
%     run(sprintf(path_EstimationCode + 'BE_EstimationCode.m'));
%     load(fullfile(path_EstimationCode, 'EstimationParameters.mat'));
% else
%     load(fullfile(path_EstimationCode, 'EstimationParameters.mat'));
% end

if exist(fullfile(path_SetupCode, 'Initial_Parameters.mat')) == 0
    run(fullfile(path_SetupCode, 'BE_DataSetup.m'));
    %load(fullfile(path_SetupCode, 'Initial_Parameters.mat'));
else
    load(fullfile(path_SetupCode, 'Initial_Parameters.mat'));
end

% da = 0.9999*(1/1.05)^(1/365);
% file = sprintf('ParameterSet%d.mat', da);
% Parameters_Table = readtable('Parameters_da0.xlsx');
% Param_Values = table2array(Parameters_Table(4,:));
% Param_Values = [12, 9, 0.0154, 0.2202, (100/(pop*1e6)), 0.1159, 0.5*0.2596, 0.5*0.2596, (1/1.05)^(1/365), 0.999*(1/1.05)^(1/365), 180, 0.1*180, 1e-12, 0.3, 0.0034, 0.3223]
Param_Values = [11, 9, 0.0154, 0.2022, (100/(pop*1e6)),  0.1159, 0.5*0.2956, 0.5*0.2956, (1/1.05)^(1/365), 180, 0.1*180, 1e-12, 0.3, 0.0034, 0.3223];
xi=20.0;
variables = {'ti', 'th', 'mi', 'mh', 'e0', 'g', 'bw', 'bs', 'dp', 'pp', 'pm', 'kap', 'Lmin', 'ifr', 'c'};
for i =1:length(variables)
    eval([variables{i} '= Param_Values(1,i);'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIALIZE VARIABLES
%%%%%%%%%%%%%%%%%%%%%
NoOfIterations = 10000;
Tolerance = 1;
Value_Function = zeros(NoOfIterations,1,length(ifr));

chi = zeros(T,1,length(ifr)); 
eta = zeros(T,T,length(ifr)); 
L = Lmin*ones(T,1,length(ifr)); 
alpha = 0.5*ones(T,T,length(ifr));  %alpha(k,t) -> k: last testing time;  t: calendar time
% load('Parameters_da99pc.mat')

SA_Final = zeros(T+1,1,length(ifr));
SR_Final = zeros(T,1,length(ifr));
IA_Final = zeros(T+1,1,length(ifr));
IR_Final = zeros(T,1,length(ifr));
tau_Final = zeros(T,1,length(ifr));
D_Final = zeros(T+1,1,length(ifr));
H_Final = zeros(T+1,1,length(ifr));

%CONVERGENCE OF VALUE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTES: 1. Tolerance level: 1
%       2. No. of iterations: 10000
%       3. Distance function: d(x,y) = sum_{i=1}^{n} |x_i-y_i|
%       4. Weighting function: 
%               (a) Iterations 0 to 9   -> 80% weight on old params; 
%               (b) Iterations 10 to 49 -> 90% weight on old params; 
%               (c) Iterations 50 to 99 -> 99% weight on old params;
%               (d) Iterations 100 to 2000 -> 99.9% weight on old params; 
tic
b=[0.99, 0.995, 0.999, 0.9995, 0.9999];
da = da.*b;
for j = 1:length(da)
da = da(j);
if exist('ParameterSet.mat')~= 0
    load('ParameterSet.mat')
end
for i=1:length(ifr)
    Figure_VF = animatedline;
    for k=1:NoOfIterations
        [chi2,eta2,alpha2,L2,SA,SR,IA,IR,tau,D, S, I, H] = fwsolve_Aug26(chi(:,:,i),eta(:,:,i),alpha(:,:,i),L(:,:,i), xi); 
        %b = 0.9*(k<10)+0.95*(k>=10)*(k<50)+0.99*(k>=50)*(k<100)+0.999*(k>=100); 
        b = b(j);
        diff_adj = sum(sum(abs(alpha(:,:,i)-alpha2)))/T+sum(abs(L(:,:,i)-L2)); 
        chi(:,:,i) = b*chi(:,:,i)+(1-b)*chi2;
        eta(:,:,i) = b*eta(:,:,i)+(1-b)*eta2;
        L(:,:,i) = b*L(:,:,i)+(1-b)*L2;
        alpha(:,:,i) = b*alpha(:,:,i)+(1-b)*alpha2;
        Value_Function(k,:,i) = diff_adj; 
        diff_adj
        if diff_adj < Tolerance
            break
        end
        addpoints(Figure_VF,k,diff_adj)
        drawnow
    end
    chi(:,:,i) = chi2;
    eta(:,:,i) = eta2;
    alpha(:,:,i) = alpha2;
    L(:,:,i) = L2;
    SA_Final(:,:,i) = SA;
    SR_Final(:,:,i) = SR;
    IA_Final(:,:,i) = IA;
    IR_Final(:,:,i) = IR;
    tau_Final(:,:,i) = tau;
    D_Final(:,:,i) = D;
    H_Final(:,:,i) = H;
end
filename = sprintf('ParameterSet%d.mat', round(da,4)*1000);
save(filename)
save('ParameterSet.mat')
end
toc

filename = sprintf('ParameterSet%d.mat', round(da,4)*1000);
save(filename)
save('ParameterSet.mat')

%SIMULATION - PLOT VALUE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OBJECTIVE: Checks if function has converged sufficiently; also checks if cycles exist
colors = {'k', 'r'};
for i=1:length(ifr)
    Figure_VF=plot(Value_Function(:,:,i), 'color', colors{i}, 'Linewidth', 1);
    hold on
    legendInfo{i} = ['ifr = ' num2str(ifr(i))];
end
title('Distance Function - Convergence (xi 20, Lmin 0.3)', 'FontSize', 13)
xlabel('Time')
ylabel('Distance between parameters')
legend(legendInfo, 'Location', 'northeast', 'FontSize' ,11)
hold off
saveas(Figure_VF, 'Value_Function_Plots/ValueFunction_xi0Lmin3.jpg')
print('Value_Function_Plots/ValueFunction_ifr34xi0Lmin3.pdf', '-dpdf', '-fillpage')

%SIMULATION - OPTIMAL CHOICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(ifr)
    Figure_OP = plot(SR_Final(:,:,i)./SA_Final(1:T,:,i),'r','LineWidth',2);
    hold on
    plot(L2(:,:,i),'k','LineWidth',2);
    title('Optimal Policy - Lockdown, Social Distancing (xi 20, Lmin 0.3)', 'FontSize', 13)
    xlabel('Time')
    ylabel('Rate of Lockdown, Rate of Contact')
    ylim([0,1.05])
    legend({ 'Eqm social distancing' 'Optimal lockdown' }, 'Location', 'southeast', 'FontSize' ,11)
    hold off
    filename = sprintf('Optimal_Policy_Plots/OptimalPolicy_ifr%dxi0Lmin3', ifr(i)*10000);
    saveas(Figure_OP, filename, 'jpg')
    print(filename, '-dpdf', '-fillpage')
end

%SIMULATION - DEATH RATES + PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Daily_Death = zeros(T,1,length(ifr));
Expected_Death = zeros(length(ifr),1);

Cumulative_Death = D_Final.*(pop*1e6);
for i=1:length(ifr)
    Daily_Death(:,:,i) = [0;diff(D_Final(1:end-1,:,i))].*(pop*1e6);
    Expected_Death(i,:) = sum(D_Final(2:end,:,i).*p)*(pop*1e6);
end

for i=1:length(ifr)
    Figure_DailyDeath = plot(Daily_Death(1:Tmax,1,i), 'r', 'Linewidth', 2);
    hold on
    plot(DNr, '--k', 'Linewidth', 1);
    title('Daily Deaths - Data and Fit (ifr 0.34%, xi 20, Lmin 0.3)', 'FontSize', 13)
    xlabel('Time')
    ylabel('No of Daily Deaths')
    legend({ 'Model' 'Data' }, 'Location', 'northwest', 'FontSize' ,11)
    hold off
    filename_DD = sprintf('Death_Rates_Plots/DailyDeaths_ifr%dxi0Lmin3', ifr(i)*10000);
    saveas(Figure_DailyDeath, filename_DD,'jpg')
    print(filename_DD, '-dpdf', '-fillpage')
end


for i=1:length(ifr)
    Figure_DailyDeath = plot(Cumulative_Death(1:T,1,i), 'r', 'Linewidth', 2);
    hold on
    plot(DCUMr, '--k', 'Linewidth', 1);
    title('Cumulative Deaths - Data and Fit (ifr 0.34%, xi 20, Lmin 0.3)', 'FontSize', 13)
    xlabel('Time')
    ylabel('Cumulative Deaths')
    legend({ 'Model' 'Data' }, 'Location', 'northwest', 'FontSize' ,11)
    hold off
    filename_DD = sprintf('Death_Rates_Plots/CumDeathsDataVModel_ifr%dxi0Lmin3', ifr(i)*10000);
    saveas(Figure_DailyDeath, filename_DD,'jpg')
    print(filename_DD, '-dpdf', '-fillpage')
end


for i=1:length(ifr)
    Figure_CumDeath = plot(Cumulative_Death(1:Tmax,1,i), 'r', 'Linewidth', 2);
    hold on
    legendInfo{i} = ['ifr = ' num2str(ifr(i))];
end
%plot(DCUMr, '--k', 'Linewidth', 1);
title('Cumulative Deaths - Data and Fit (xi 20, Lmin 0.3)', 'FontSize', 13)
xlabel('Time')
ylabel('No of Cumulative Deaths')
legend(legendInfo, 'Location', 'northwest', 'FontSize' ,11)
hold off
filename_CD = sprintf('Death_Rates_Plots/CumDeaths_ifrxi0Lmin3');
saveas(Figure_CumDeath, filename_CD, 'jpg')
print(filename_CD, '-dpdf', '-fillpage')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













