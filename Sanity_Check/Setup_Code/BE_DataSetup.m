%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                           BEHAVIOURAL EPIDEMIOLOGY
%                                 SETUP CODE
%                                SANITY CHECK
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
global da T Tmax p q X Lmin pop POSN DN HST

%PATHNAMES REQUIRED - change as required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_SetupPLots = '../Figures/Setup_Plots/';  %This is the folder where the miscellaneous plots are saved - vaccine arrival probability, daily tests fit
path_SetupParameters  = '../Parameters/BE_DataSetup/';
path_OldParameters = '../Parameters/Old_Parameters/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEFINE PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%
Tm = 540;
Tv = 180;
pop = 328.2;
p = zeros(T,1);
q =  zeros(T,1);
cdf_p = zeros(T,1);

%Create a negative binomial distribution: parameters -> r,p. mean = r/p; variance = [r(1-p)]/p^2
Ns = ceil(Tm*Tm/(Tm+Tv));
f = makedist('NegativeBinomial','R',Ns,'p',Tm/(Tm+Tv));
Nt = icdf(f,1-1e-4);
T = Ns+Nt;
p(Ns,1) = cdf(f,0);
for n=Ns+1:T
    p(n) = cdf(f,n-Ns+1)-cdf(f,n-Ns);
end
for n=1:T
    q(n) = sum(p(n:T));
end
for n=1:T
    cdf_p(n) = sum(p(1:n));
end

%DEFINE PROBABILITIES - MISCELLANEOUS PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figure_p = plot(p,'^-k', 'Linewidth', 1);
title('Probability of vaccine arrival', 'FontSize', 13)
xlabel('Time')
ylim([0,0.04])
ylabel('Probability')
legend({ 'p(t) = prob. of arrival at t' }, 'Location', 'northwest', 'FontSize' ,12)
saveas(Figure_p, fullfile(path_SetupPLots, 'pdf_VaccineArrival'), 'jpeg');
saveas(Figure_p, fullfile(path_SetupPLots, 'pdf_VaccineArrival'), 'pdf');

Figure_cdf_p = plot(cdf_p, 'd-k', 'Linewidth', 1);
title('Cumulative probability of vaccine arrival', 'FontSize', 13)
xlabel('Time')
ylim([0,1.05])
ylabel('Probability')
legend({ 'cdf_p(t) = Pr[vaccine arrival time<=t]' }, 'Location', 'northwest', 'FontSize' ,12)
saveas(Figure_cdf_p, fullfile(path_SetupPLots, 'cdf_VaccineArrival'), 'jpeg');
saveas(Figure_cdf_p, fullfile(path_SetupPLots, 'cdf_VaccineArrival'), 'pdf');

Figure_q = plot(q, 's-k', 'Linewidth', 1);
title('Probability of no vaccine', 'FontSize', 13)
xlabel('Time')
ylim([0,1.05])
ylabel('Probability')
legend({ 'q(t) = Pr[vaccine arrival time>=t]' }, 'Location', 'southwest', 'FontSize' ,12)
saveas(Figure_q, fullfile(path_SetupPLots, 'q_NoVaccineArrival'), 'jpeg');
saveas(Figure_q, fullfile(path_SetupPLots, 'q_NoVaccineArrival'), 'pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%READ DATASET
%%%%%%%%%%%%%
%The dataset will only be required to fit tests X(t) and plot graphs - simulations do not require the data
data = readtable('../Dataset/data_upd3xlsx.xlsx');

%READ VARIABLES
%%%%%%%%%%%%%%%
POSCUMr = data.POSCUM;
NEGCUMr = data.NEGCUM;
DCUMr = data.DCUM;
RCUMr = data.RCUM;
Hr = data.H;

%replace all missing values (recorded as NaN's)  with 0's
POSCUMr(isnan(POSCUMr))=0;
NEGCUMr(isnan(NEGCUMr))=0;
DCUMr(isnan(DCUMr))=0;
RCUMr(isnan(RCUMr))=0;
Hr(isnan(Hr)) = 0;

%NOTE: Shift each variable one day forward, because the data is recorded at the end of each day
POSNr = [0;diff(POSCUMr(1:end-1))]; %# of NEW positive tests at t (raw)
NEGNr = [0;diff(NEGCUMr(1:end-1))]; %# of NEW negative tests at t (raw)
DNr = [0;diff(DCUMr(1:end-1))]; %# of NEW fatal cases at t (raw)
RNr = [0; diff(RCUMr(1:end-1))]; %# of NEW recovered cases at t (raw)

%smoothing the functions out
POSN = smooth(POSNr,5,'moving');
NEGN = smooth(NEGNr,5,'moving');
DN = smooth(DNr,5,'moving');
RN = smooth(RNr,5,'moving');
HST = smooth(Hr,5,'moving');
HST = HST(1:end-1);

%Calculate Optimal Weighting Function
Weight = inv(cov([POSNr NEGNr DNr Hr(1:end-1)]));

%CURVE FIT TEST VARIABLE FROM DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: Total tests at time t = No. of positive tests at time t + No. of negative tests at time t

X = ones(T,1)/(pop*1e6);
Tmax = 131 - 1;
%Tmax = length(data.H)-1;

%record current calendar time, that is x(t) = t
x = 1:1:T;
x = x';

%Linear fit for tests between time Tr to Tm: (POSNr+NEGNr)(t) = beta_0+beta_1*x(t)
fitx = fit(x(1:Tmax-1),POSNr(1:Tmax-1)+NEGNr(1:Tmax-1),'poly2');
for t=1:T
    X(t) = 1.*min(1e6,1*max(1,fitx(x(t))))/(pop*1e6);
end

%PLOT X(t) AND TEST DATA
%%%%%%%%%%%%%%%%%%%%%%%%
Figure_Test = plot(X.*(pop*1e6), 'r', 'Linewidth', 2);
hold on
plot(POSNr+NEGNr, '--k', 'Linewidth', 1);
title('Daily Tests - Data and Fit', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Tests')
ylim([0,11*1e5])
legend({ 'Model' 'Data' }, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_Test, fullfile(path_SetupPLots, 'DailyTests_Data'), 'jpeg');
saveas(Figure_Test, fullfile(path_SetupPLots, 'DailyTests_Data'), 'pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT DATA
%%%%%%%%%%
Figure_DNr = plot(DNr(1:Tmax), '--k', 'Linewidth', 1);
title('Daily Deaths - Data', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Deaths')
legend({ 'Data'}, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_DNr, fullfile(path_SetupPLots, 'DailyDeaths_Data'), 'jpeg');
saveas(Figure_DNr, fullfile(path_SetupPLots, 'DailyDeaths_Data'), 'pdf');

Figure_HST = plot(HST(1:Tmax), '--k', 'Linewidth', 1);
title('Daily Hospitalized - Data', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Hospitalized')
legend({ 'Data'}, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_HST, fullfile(path_SetupPLots, 'DailyHospitalized_Data'), 'jpeg');
saveas(Figure_HST, fullfile(path_SetupPLots, 'DailyHospitalized_Data'), 'pdf');

Figure_NEGNr = plot(NEGNr(1:Tmax), '--k', 'Linewidth', 1);
title('Daily Negative Tests - Data', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Negative Tests')
legend({ 'Data'}, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_NEGNr, fullfile(path_SetupPLots, 'DailyNegCases_Data'), 'jpeg');
saveas(Figure_NEGNr, fullfile(path_SetupPLots, 'DailyNegCases_Data'), 'pdf');

Figure_POSNr = plot(POSNr(1:Tmax), '--k', 'Linewidth', 1);
title('Daily Positive Cases - Data', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Positive Cases')
legend({ 'Data'}, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_POSNr, fullfile(path_SetupPLots, 'DailyPosCases_Data'), 'jpeg');
saveas(Figure_POSNr, fullfile(path_SetupPLots, 'DailyPosCases_Data'), 'pdf');

Figure_RNr = plot(RNr(1:Tmax), '--k', 'Linewidth', 1);
title('Daily Recoveries - Data', 'FontSize', 14)
xlabel('Time')
ylabel('No of Daily Recoveries')
legend({ 'Data'}, 'Location', 'southeast', 'FontSize' ,11)
hold off
saveas(Figure_RNr, fullfile(path_SetupPLots, 'DailyRecoveries_Data'), 'jpeg');
saveas(Figure_RNr, fullfile(path_SetupPLots, 'DailyRecoveries_Data'), 'pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETER VALUES
%%%%%%%%%%%%%%%%%
%NOTE: These are the fixed variables.

%da = 0.99*(1/1.05)^(1/365);
da = (1/1.05)^(1/365);
Fixed_ParamValues = [10.4685, 9.2312, (100/(pop*1e6)), (1/1.05)^(1/365), 1e-12, 0.4, 0.002, 5, 0.5];
Fixed_Variables = {'ti', 'th', 'e0', 'dp', 'kap', 'Lmin', 'ifr', 'rl1', 'beta_share'};
for i = 1:length(Fixed_Variables)
    eval([Fixed_Variables{i} '= Fixed_ParamValues(1,i);'])
end
% mi = (Hr(1:Tmax)\DNr(1:Tmax))*th;
% mh = ifr/mi;
mi = 0.0076;
mh = 0.2618;
xi=20;

%UPDATE THE PARAMETERS FROM REGRESSIONS RUN BY ILIA FOR MYOPIC AGENTS
% Parameters_Table = readtable(fullfile(path_OldParameters,'Parameters_da0.xlsx'));
% Param_Values = table2array(Parameters_Table(4,:));
% variables = {'ti', 'th', 'mi', 'mh', 'e0', 'g', 'bw', 'bs', 'dp', 'da', 'pp', 'pm', 'kap', 'Lmin', 'ifr', 'c'};
% for i =1:length(variables)
%     eval([variables{i} '= Param_Values(1,i);'])
% end
% xi=20;

filename = sprintf('Initial_Parameters');
save(fullfile(path_SetupParameters,filename));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
