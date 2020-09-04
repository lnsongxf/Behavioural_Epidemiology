function [chi2,eta2,alpha2,L2,SA,SR,IA,IR,tau,D, S, I, H] = fwsolve_Aug26(chi,eta,alpha,L, xi)

%pp = phi_plus// instantenous costs of infecting another person
%pm = phi_minus// instantenous costs of getting infected
%c = flow costs of social distancing
%g = gamma// efficiency tracing
%da = delta_a// agents' discount factor
%dp = delta_p// government's discount factor
%T = number of periods the model is solved for
%p = pdf of vaccine arrival time, i.e. p(t) = probability it arrives at t
%q = probability that the vaccine has not arrived until t, that is q(t) =sum_{k >= t}p(k)
%e0 = fraction of people infected at the outset
%bw = beta_w// contact rate in productive activities
%bs = beta_s// contact rate in social activities
%ti = average length of a mild infection
%th = average length of hospital stay
%mi = hospitalization rate
%mh = death rate
%X = testing capacity, i.e. X(t) = number of tests available at t
%Lmin = lower bound on lockdown
%kap = kappa, tunning parameter to smooth the agent's social distancing problem

global pp pm c g da dp T p q e0 bw bs ti th mi mh X Lmin kap

%%%%%% SYSTEM I.
%%%%%%%%%%%%%%%%

%initialize all state variables where
%t - current calendar time, k - last testing date
S = zeros(T+1,T+1); %stock of suscetible, S(k,t), matrix of size T+1, T+1
I = zeros(T+1,T+1); %stock of infected (unknown), S(k,t), matrix of size T+1, T+1
R = zeros(T+1,T+1); %stock of recovered (unknown), S(k,t), matrix of size T+1, T+1
SR = zeros(T,1); %stock of susceptible participating in in social activities, SR(t) = sum_{k<=t}alpha_t^kS_t^k, vector of length T
IR = zeros(T,1); %stock of infected (unknown) participating in in social activities, IR(t) = sum_{k<=t}alpha_t^kI_t^k, vector of length T
SA = zeros(T+1,1); %total stock of susceptible, SA(t) = sum_{k<=t}S_t^k, vector of length T+1
IA = zeros(T+1,1); %total stock of infected (unknown), IA(t) = sum_{k<=t}I_t^k, vector of length T+1
RA = zeros(T+1,1); %total stock of recovered (unknown), RA(t) = sum_{k<=t}R_t^k, vector of length T+1
IT = zeros(T+1,1); %stock of infected (known), vector of length T+1
RT = zeros(T+1,1); %stock of recovered (known), vector of length T+1
H = zeros(T+1,1); %stock of hospitalized, vector of length T+1
D = zeros(T+1,1); %stock of dead, vector of length T+1
tau = zeros(T,1); %testing probability, vector of length T
S(1,1) = 1-e0; %fraction of susceptible at the outset is 1-e0
SA(1) = 1-e0; 
I(1,1) = e0; %fraction of infected (known) at the outset is e0
IA(1) = e0;

%given alpha,L sovle the system of state forward equations
for t=1:T    
    SR(t) = sum(alpha(1:t,t).*S(1:t,t)); %find SR(t) = sum_{k<=t}alpha_t^kS_t^k
    IR(t) = sum(alpha(1:t,t).*I(1:t,t)); %find IR(t) = sum_{k<=t}alpha_t^kI_t^k
    SA(t+1) = SA(t)-bw*(L(t))^2*SA(t)*IA(t)-bs*SR(t)*IR(t); 
    tau(t) = X(t)/(g*SA(t+1)+g*(RA(t)+(1-mi)/ti*IA(t))+(1-1/ti)*IA(t)+bw*(L(t))^2*SA(t)*IA(t)+bs*SR(t)*IR(t)); 
    IA(t+1) = (1-tau(t))*((1-1/ti)*IA(t)+bw*(L(t))^2*SA(t)*IA(t)+bs*SR(t)*IR(t));
    RA(t+1) = (1-tau(t)*g)*(RA(t)+(1-mi)/ti*IA(t)); 
    S(t+1,t+1) = tau(t)*g*SA(t+1);
    IT(t+1) = (1-1/ti)*IT(t)+tau(t)*((1-1/ti)*IA(t)+bw*(L(t))^2*SA(t)*IA(t)+bs*SR(t)*IR(t)); 
    RT(t+1) = RT(t)+(1-mi)/ti*IT(t)+tau(t)*g*(RA(t)+(1-mi)/ti*IA(t))+(1-mh)/th*H(t); 
    H(t+1) = (1-1/th)*H(t)+mi/ti*(IA(t)+IT(t)); 
    %update S(k,t+1),I(k,t+1),R(k,t+1) for k=1,...,t
    S(1:t,t+1)=(1-g*tau(t))*(S(1:t,t)-bw*(L(t))^2*S(1:t,t)*IA(t)-bs*alpha(1:t,t).*S(1:t,t)*IR(t));
    I(1:t,t+1)=(1-tau(t))*(I(1:t,t)*(1-1/ti)+bw*(L(t))^2*S(1:t,t)*IA(t)+bs*alpha(1:t,t).*S(1:t,t)*IR(t));    
    R(1:t,t+1)=(1-g*tau(t))*(R(1:t,t)+(1-mi)/ti*I(1:t,t));
    D(t+1) = D(t) + mh/th*H(t); %compute D(t+1)
end

%%%%%% SYSTEM II.
%%%%%%%%%%%%%%%%

%initialize all adjoint state variables for the agents' adjoint equations, here 
%the interpretations as before and '2' means that the variable is the adjoint to the agents' adjoints

S2 = zeros(T+1,T+1); 
I2 = zeros(T+1,T+1);
R2 = zeros(T+1,T+1);
SR2 = zeros(T,1);
IR2 = zeros(T,1);
SA2 = zeros(T+1,1);
IA2 = zeros(T+1,1);
RA2 = zeros(T+1,1);
IT2 = zeros(T+1,1);
RT2 = zeros(T+1,1);
H2 = zeros(T+1,1);

%initial conditions at t=2
S2(1,2) = -(1-tau(1)*g)*eta(1,1)*bs*SA(1)*IR(1);
S2(2,2) = -tau(1)*g*eta(1,1)*bs*SA(1)*IR(1);
SA2(2) = S2(1,2) + S2(2,2);
I2(1,2) = (1-tau(1))*eta(1,1)*bs*SA(1)*IR(1);
IA2(2) = I2(1,2);
IT2(2) = tau(1)*eta(1,1)*bs*SA(1)*IR(1);

%given alpha,L,eta sovle the system of agents' forward equations
for t=2:T-1    
    SR2(t) = sum(alpha(1:t,t).*S2(1:t,t));  %find SR2(t) = sum_{k<=t}alpha_t^k \bar {S}_t^k
    IR2(t) = sum(alpha(1:t,t).*I2(1:t,t)); %find IR2(t) = sum_{k<=t}alpha_t^k \bar {I}_t^k
    SA2(t+1) = da*(SA2(t)-bw*(L(t))^2*SA2(t)*IA(t)-bs*SR2(t)*IR(t))...
        -bs*sum(eta(1:t,t).*S(1:t,t))*IR(t); 
    IA2(t+1) = da*(1-tau(t))*((1-1/ti)*IA2(t)+bw*(L(t))^2*SA2(t)*IA(t)+bs*SR2(t)*IR(t))...
        +(1-tau(t))*bs*sum(eta(1:t,t).*S(1:t,t))*IR(t); 
    RA2(t+1) = da*(1-tau(t)*g)*(RA2(t)+(1-mi)/ti*IA2(t)); 
    S2(t+1,t+1) = tau(t)*g*SA2(t+1);
    IT2(t+1) = da*((1-1/ti)*IT2(t)+tau(t)*((1-1/ti)*IA2(t)+bw*(L(t))^2*SA2(t)*IA(t)+bs*SR2(t)*IR(t)))...
        +tau(t)*bs*sum(eta(1:t,t).*S(1:t,t))*IR(t);
    RT2(t+1) = da*(RT2(t)+(1-mi)/ti*IT2(t)+tau(t)*g*(RA2(t)+(1-mi)/ti*IA2(t))+(1-mh)/th*H2(t));
    H2(t+1) = da*((1-1/th)*H2(t)+mi/ti*(IA2(t)+IT2(t)));
    %update S2(k,t+1),I2(k,t+1),R2(k,t+1) for k=1,...,t    
    S2(1:t,t+1)= da*(1-g*tau(t))*(S2(1:t,t)-bw*(L(t))^2*S2(1:t,t)*IA(t)-bs*alpha(1:t,t).*S2(1:t,t)*IR(t))...
        -(1-tau(t)*g)*bs*eta(1:t,t).*S(1:t,t)*IR(t);
    I2(1:t,t+1)= da*(1-tau(t))*(I2(1:t,t)*(1-1/ti)+bw*(L(t))^2*S2(1:t,t)*IA(t)+bs*alpha(1:t,t).*S2(1:t,t)*IR(t))...
        +(1-tau(t))*bs*eta(1:t,t).*S(1:t,t)*IR(t);
    R2(1:t,t+1)= da*(1-tau(t)*g)*(R2(1:t,t)+(1-mi)/ti*I2(1:t,t));  
end 

%%%%%% SYSTEM III.
%%%%%%%%%%%%%%%%

    
%initialize all agents' adjoint variables where
%t - current calendar time, k - last testing date
Sx = zeros(T+1,T); %multiplier on S(k,t+1) = ..., T+1,T matrix
Sx(:,T) = da*(c/2+1)*p(T)*ones(T+1,1); %terminal conditions for Sx
Ix = zeros(T+1,T); %multiplier on I(k,t+1) = ..., T+1,T matrix
Ix(1:T,T) = da*(c/2+1)*p(T)*ones(T,1); %terminal conditions for Ix
ITx = zeros(1,T); %multiplier on IT(t+1) = ..., 1,T matrix
ITx(:,T) = da*(c/2+1)*p(T); %terminal conditions for ITx
Rx = zeros(T+1,T); %multiplier on R(k,t+1) = ..., T+1,T matrix
Rx(1:T,T) = da*(c/2+1)*p(T)*ones(T,1); %terminal conditions for Rx
RTx = zeros(1,T); %multiplier on RT(t+1) = ..., 1,T matrix
RTx(:,T) = da*(c/2+1)*p(T); %terminal conditions for RTx
Hx = zeros(1,T); %multiplier on H(t+1) = ..., 1,T matrix
Hx(:,T) = da*(c/2+1)*p(T); %terminal conditions for Hx

%given the forward variables,tau,alpha,L solve the agents' backward system, here
%t=1 is one day before the terminal date and t=T-1 means the first date
for t=1:T-1
    Sx(1:T+1-t,T-t) = da*((1-g*tau(T+1-t))*Sx(1:T+1-t,T+1-t)+g*tau(T+1-t)*Sx(T+2-t,T+1-t))...
        +da*((1-tau(T+1-t))*Ix(1:T+1-t,T+1-t)+tau(T+1-t)*ITx(T+1-t)-(1-g*tau(T+1-t))*Sx(1:T+1-t,T+1-t)-g*tau(T+1-t)*Sx(T+2-t,T+1-t)).*(bw*(L(T+1-t))^2*IA(T+1-t)+bs*alpha(1:T+1-t,T+1-t)*IR(T+1-t))...
        +da*(1-da)*L(T+1-t)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha(1:T+1-t,T+1-t)).^2*(c/2)*q(T+1-t)...
        -da*(1-da)*(L(T+1-t))^2*IA(T+1-t).*bw.*pp*q(T+1-t)...
        -da*(1-da)*alpha(1:T+1-t,T+1-t)*IR(T+1-t).*bs.*pp*q(T+1-t);
    Ix(1:T-t,T-t) = da*(1-1/ti)*((1-tau(T+1-t))*Ix(1:T-t,T+1-t)+tau(T+1-t)*ITx(T+1-t))...
        +da*(1-mi)/ti*((1-g*tau(T+1-t))*Rx(1:T-t,T+1-t)+g*tau(T+1-t)*RTx(T+1-t))+da*mi/ti*Hx(T+1-t)...
        +da*(1-da)*L(T+1-t)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha(1:T-t,T+1-t)).^2*(c/2)*q(T+1-t)...
        -da*(1-da)*(L(T+1-t))^2*SA(T+1-t).*bw*pm*q(T+1-t)...
        -da*(1-da)*alpha(1:T-t,T+1-t)*SR(T+1-t).*bs*pm*q(T+1-t);
    Rx(1:T-t,T-t) = da*((1-g*tau(T+1-t))*Rx(1:T-t,T+1-t)+g*tau(T+1-t)*RTx(T+1-t))...
        +da*(1-da)*L(T+1-t)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha(1:T-t,T+1-t)).^2*(c/2)*q(T+1-t);
    ITx(T-t) = da*((1-1/ti)*ITx(T+1-t)+(1-mi)/ti*RTx(T+1-t)+mi/ti*Hx(T+1-t))+da*p(T-t)+da*(c/2)*p(T-t);
    RTx(T-t) = da*RTx(T+1-t)+da*(1-da)*L(T+1-t)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t);
    Hx(T-t) = da*((1-1/th)*Hx(T+1-t)+(1-mh)/th*RTx(T+1-t))+da*p(T-t)+da*(c/2)*p(T-t);
end
   
%%%%%% SYSTEM IV.
%%%%%%%%%%%%%%%%

%initialize all government' adjoint variables here 
%the interpretations as before and '2' means that the variable is the government's adjoint
Sx2 = zeros(T+1,T);
Sx2(:,T) = (xi+dp^T)*p(T)*ones(T+1,1);
Ix2 = zeros(T+1,T);
Ix2(1:T,T) = (xi+dp^T)*p(T)*ones(T,1);
ITx2 = zeros(1,T);
ITx2(:,T) = (xi+dp^T)*p(T);
Rx2 = zeros(T+1,T);
Rx2(1:T,T) = (xi+dp^T)*p(T)*ones(T,1);
RTx2 = zeros(1,T);
RTx2(:,T) = (xi+dp^T)*p(T);
Hx2 = zeros(1,T);
Hx2(:,T) = (xi+dp^T)*p(T);

for t=1:T-1
    
    Sx2(1:T+1-t,T-t)  =  ((1-g*tau(T+1-t))*Sx2(1:T+1-t,T+1-t)+g*tau(T+1-t)*Sx2(T+2-t,T+1-t))...
         +((1-tau(T+1-t))*Ix2(1:T+1-t,T+1-t)+tau(T+1-t)*ITx2(T+1-t)-(1-g*tau(T+1-t))*Sx2(1:T+1-t,T+1-t)-g*tau(T+1-t)*Sx2(T+2-t,T+1-t)).*(bw*(L(T+1-t))^2*IA(T+1-t)+bs*alpha(1:T+1-t,1)*IR(T+1-t))...
         +eta(1:T+1-t,T+1-t).*(c*(1-alpha(1:T+1-t,T+1-t))-bs.*pp*IR(T+1-t))*(1-da)*q(T+1-t)-bs.*pm*sum(eta(1:T+1-t,T+1-t).*I(1:T+1-t,T+1-t))*alpha(1:T+1-t,T+1-t)*(1-da)*q(T+1-t)...
         +((1-tau(T+1-t))*Ix(1:T+1-t,T+1-t)+tau(T+1-t)*ITx(T+1-t)-(1-tau(T+1-t)*g)*Sx(1:T+1-t,T+1-t)-tau(T+1-t)*g*Sx(T+2-t,T+1-t))*bs.*eta(1:T+1-t,T+1-t)*IR(T+1-t)...
         -chi(T+1-t)*tau(T+1-t)*(g+(1-g)*(bw*(L(T+1-t))^2*IA(T+1-t)+bs*alpha(1:T+1-t,T+1-t)*IR(T+1-t)))...
         -da*(1-da)*IA2(T+1-t)*(L(T+1-t))^2.*bw*pm*q(T+1-t)...
         -da*(1-da)*IR2(T+1-t)*alpha(1:T+1-t,T+1-t).*bs*pm*q(T+1-t)...                  
         +(1-dp)*(dp)^(T-t)*L(T+1-t)*q(T+1-t)+(xi+dp^(T-t))*p(T-t);
     
    Ix2(1:T-t,T-t)  = (1-1/ti)*((1-tau(T+1-t))*Ix2(1:T-t,T+1-t)+tau(T+1-t)*ITx2(T+1-t))...
        +(1-mi)/ti*((1-g*tau(T+1-t))*Rx2(1:T-t,T+1-t)+g*tau(T+1-t)*RTx2(T+1-t))+mi/ti*Hx2(T+1-t)...
        +sum(((1-tau(T+1-t))*Ix2(1:T+1-t,T+1-t)-(1-tau(T+1-t)*g)*Sx2(1:T+1-t,T+1-t)).*S(1:T+1-t,T+1-t))*bw*(L(T+1-t))^2 ...
        +sum(((1-tau(T+1-t))*Ix2(1:T+1-t,T+1-t)-(1-tau(T+1-t)*g)*Sx2(1:T+1-t,T+1-t)).*alpha(1:T+1-t,T+1-t).*S(1:T+1-t,T+1-t))*bs.*alpha(1:T-t,T+1-t)...
        +(tau(T+1-t)*ITx2(T+1-t)-tau(T+1-t)*g*Sx2(T+2-t,T+1-t))*(bw*(L(T+1-t))^2*SA(T+1-t)+bs*SR(T+1-t).*alpha(1:T-t,T+1-t))...
        +da*sum(((1-tau(T+1-t))*Ix(1:T+1-t,T+1-t)-(1-tau(T+1-t)*g)*Sx(1:T+1-t,T+1-t)).*S2(1:T+1-t,T+1-t))*bw*(L(T+1-t))^2 ...
        +da*sum(((1-tau(T+1-t))*Ix(1:T+1-t,T+1-t)-(1-tau(T+1-t)*g)*Sx(1:T+1-t,T+1-t)).*alpha(1:T+1-t,T+1-t).*S2(1:T+1-t,T+1-t))*bs.*alpha(1:T-t,T+1-t)...
        +da*(tau(T+1-t)*ITx(T+1-t)-tau(T+1-t)*g*Sx(T+2-t,T+1-t))*(bw*(L(T+1-t))^2*SA2(T+1-t)+bs*SR2(T+1-t).*alpha(1:T-t,T+1-t))...
        -da*(1-da)*SA2(T+1-t)*(L(T+1-t))^2*bw*pp*q(T+1-t)...
        -da*(1-da)*SR2(T+1-t)*alpha(1:T-t,T+1-t)*bs*pp*q(T+1-t)...
        +(1-dp)*(dp)^(T-t)*L(T+1-t)*q(T+1-t)+(xi+dp^(T-t))*p(T-t)...
        +eta(1:T-t,T+1-t).*(c*(1-alpha(1:T-t,T+1-t))-bs.*pm*SR(T+1-t))*(1-da)*q(T+1-t)-bs.*pp*sum(eta(1:T+1-t,T+1-t).*S(1:T+1-t,T+1-t))*alpha(T-t,T+1-t)*(1-da)*q(T+1-t)...
        +sum(((1-tau(T+1-t))*Ix(1:T+1-t,T+1-t)+tau(T+1-t)*ITx(T+1-t)-(1-tau(T+1-t)*g)*Sx(1:T+1-t,T+1-t)-tau(T+1-t)*g*Sx(T+2-t,T+1-t)).*S(1:T+1-t,T+1-t).*eta(1:T+1-t,T+1-t))*bs*alpha(1:T-t,T+1-t)...
        -chi(T+1-t)*tau(T+1-t)*(1-1/ti+g*(1-mi)/ti+(1-g)*(bw*(L(T+1-t))^2*SA(T+1-t)+bs*alpha(1:T-t,T+1-t)*SR(T+1-t)));
    
    Rx2(1:T-t,T-t) = ((1-g*tau(T+1-t))*Rx2(1:T-t,T+1-t)+g*tau(T+1-t)*RTx2(T+1-t))...
        +eta(1:T-t,T+1-t)*c.*(1-alpha(1:T-t,T+1-t))*(1-da)*q(T+1-t)-chi(T+1-t)*tau(T+1-t)*g+(1-dp)*(dp)^(T-t)*L(T+1-t)*q(T+1-t)+(xi+dp^(T-t))*p(T-t);
    
    ITx2(T-t) = (1-1/ti)*ITx2(T+1-t)+(1-mi)/ti*RTx2(T+1-t)+mi/ti*Hx2(T+1-t)+(xi+dp^(T-t))*p(T-t);
    RTx2(T-t) = RTx2(T+1-t)+(1-dp)*(dp)^(T-t)*L(T+1-t)*q(T+1-t)+(xi+dp^(T-t))*p(T-t);
    Hx2(T-t) = ((1-1/th)*Hx2(T+1-t)+(1-mh)/th*RTx2(T+1-t))+(xi+dp^(T-t))*p(T-t);
end
    
%initialize updated policies
chi2 = zeros(T,1); %multiplier on testing,  T vector
eta2 = zeros(T,T); %multiplier on social distancing foc,  T,T matrix
L2 = ones(T,1)*1; %lockdown,  T vector
alpha2 = zeros(T,T); %social distancing T,T matrix


for t=1:T
    %update chi // see the note for this formula
    chi2(t)= -g*sum(Sx2(1:t,t).*(S(1:t,t)-bw*(L(t))^2*S(1:t,t)*IA(t)-bs*alpha(1:t,t).*S(1:t,t)*IR(t)))...
        +g*Sx2(t+1,t)*(SA(t)-bw*(L(t))^2*SA(t)*IA(t)-bs*SR(t)*IR(t))...
        -sum(Ix2(1:t,t).*(I(1:t,t)*(1-1/ti)+bw*(L(t))^2*S(1:t,t)*IA(t)+bs*alpha(1:t,t).*S(1:t,t)*IR(t)))...
        +ITx2(t)*((1-1/ti)*IA(t)+bw*(L(t))^2*SA(t)*IA(t)+bs*SR(t)*IR(t))...
        -g*sum(Rx2(1:t,t).*(R(1:t,t)+(1-mi)/ti*I(1:t,t)))...
        +g*RTx2(t)*(RA(t)+(1-mi)/ti*IA(t))...
        -da*g*sum(Sx(1:t,t).*(S2(1:t,t)-bw*(L(t))^2*S2(1:t,t)*IA(t)-bs*alpha(1:t,t).*S2(1:t,t)*IR(t)))...
        +da*g*Sx(t+1,t)*(SA2(t)-bw*(L(t))^2*SA2(t)*IA(t)-bs*SR2(t)*IR(t))...
        -da*sum(Ix(1:t,t).*(I2(1:t,t)*(1-1/ti)+bw*(L(t))^2*S2(1:t,t)*IA(t)+bs*alpha(1:t,t).*S2(1:t,t)*IR(t)))...
        +da*ITx(t)*((1-1/ti)*IA2(t)+bw*(L(t))^2*SA2(t)*IA(t)+bs*SR2(t)*IR(t))...
        -da*g*sum(Rx(1:t,t).*(R2(1:t,t)+(1-mi)/ti*I2(1:t,t)))...
        +da*g*RTx(t)*(RA2(t)+(1-mi)/ti*IA2(t))...
        -sum(eta(1:t,t)*bs.*Sx(1:t,t)*IR(t).*(g*S(t+1,t)-g*Sx(1:t,t)-ITx(t)+Ix(1:t,t)));
    
   chi2(t) = chi2(t)/(g*SA(t+1)+g*(RA(t)+(1-mi)/ti*IA(t))+((1-1/ti)*IA(t)+(bw*(L(t))^2*SA(t)*IA(t)+bs*SR(t)*IR(t))));
    
  
    %update lockdown // see the note for this formula
    a = (1-dp)*(dp)^(t-1)*(SA(t)+IA(t)+RA(t)+RT(t))*q(t)+(1-da)*da*(SA2(t)+IA2(t)+RA2(t)+RT2(t))*q(t);
    b = (1-g)*chi(t)*bw*tau(t)*SA(t)*IA(t)...
        -bw*sum(((1-tau(t))*Ix2(1:t,t)+tau(t)*ITx2(t)-(1-tau(t)*g)*Sx2(1:t,t)-tau(t)*g*Sx2(t+1,t)).*S(1:t,t))*IA(t)...
        -da*bw*sum(((1-tau(t))*Ix(1:t,t)+tau(t)*ITx(t)-(1-tau(t)*g)*Sx(1:t,t)-tau(t)*g*Sx(t+1,t)).*S2(1:t,t))*IA(t)...
        +da*(1-da)*bw*pp*SA2(t)*IA(t)*q(t)+da*(1-da)*bw*pm*SA(t)*IA2(t)*q(t);
    if b<=0
        L2(t)=1;
    else
        L2(t)=min(1,max(Lmin,a/(2*b)));
    end
    
    %update eta // see the note for this formula
    A = bs*(1-da)*(pp.*S(1:t,t)'.*I(1:t,t)+pm.*I(1:t,t)'.*S(1:t,t))*q(t)...
        -bs*((1-tau(t))*Ix(1:t,t)'+tau(t)*ITx(t)-(1-tau(t)*g)*Sx(1:t,t)'-tau(t)*g*Sx(t+1,t)).*S(1:t,t)'.*I(1:t,t)...
        +(1-da)*c*diag(S(1:t,t)+I(1:t,t)+R(1:t,t)+kap*(1./(alpha(1:t,t)).^2+1./(1-alpha(1:t,t)).^2))*q(t);
    
    B = bs*((1-tau(t))*Ix2(1:t,t)+tau(t)*ITx2(t)-(1-tau(t)*g)*Sx2(1:t,t)-tau(t)*g*Sx2(t+1,t)).*S(1:t,t)*IR(t)...
        +da*bs*((1-tau(t))*Ix(1:t,t)+tau(t)*ITx(t)-(1-tau(t)*g)*Sx(1:t,t)-tau(t)*g*Sx(t+1,t)).*S2(1:t,t)*IR(t)...
        +bs*sum(((1-tau(t))*Ix2(1:t,t)+tau(t)*ITx2(t)-(1-tau(t)*g)*Sx2(1:t,t)-tau(t)*g*Sx2(t+1,t)).*alpha(1:t,t).*S(1:t,t)).*I(1:t,t)...
        +da*bs*sum(((1-tau(t))*Ix(1:t,t)+tau(t)*ITx(t)-(1-tau(t)*g)*Sx(1:t,t)-tau(t)*g*Sx(t+1,t)).*alpha(1:t,t).*S2(1:t,t)).*I(1:t,t)...        
        +da*(1-da)*c*(1-alpha(1:t,t)).*(S2(1:t,t)+I2(1:t,t)+R2(1:t,t))*q(t)...
        -da*(1-da)*(bs*pp.*(S2(1:t,t)*IR(t)+I(1:t,t)*SR2(t))+bs*pm.*(I2(1:t,t)*SR(t)+S(1:t,t)*IR2(t)))*q(t)...
        -(1-g)*bs*tau(t)*chi(t)*(S(1:t,t)*IR(t)+I(1:t,t)*SR(t));
    eta2(1:t,t) = A\B;
    
    %update alpha // see the note for this formula
    a = 1-bs*pp/c*S(1:t,t)./(1e-24+S(1:t,t)+I(1:t,t)+R(1:t,t))*IR(t)-bs*pm/c*I(1:t,t)./(1e-24+S(1:t,t)+I(1:t,t)+R(1:t,t))*SR(t)+bs/c*1/(q(t)*(1-da))*S(1:t,t)./(1e-24+S(1:t,t)+I(1:t,t)+R(1:t,t))*IR(t).*((1-tau(t))*Ix(1:t,t)+tau(t)*ITx(t)-(1-tau(t)*g)*Sx(1:t,t)-tau(t)*g*Sx(t+1,t));
%     b = kap./(1e-24+S(1:t,t)+I(1:t,t)+R(1:t,t));
%     a1 = (3*(a-2*b)-(1+a).^2)/3;
%     b1 = (9*(1+a).*(a-2*b)-2*(1+a).^3+27*b)/27;
%     alpha2(1:t,t) = (1+a)/3+2*sqrt(-a1/3).*cos(acos(3/2*b1./a1.*sqrt(-3./a1))/3-2/3*pi);
    alpha2(1:t,t) = max(0, min(1,a));
end
end



