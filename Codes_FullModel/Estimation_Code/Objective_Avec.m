function [SS, XP0, XN0, D0, SA, SR, L, tau, alpha] = Objective_Avec(pars)

global ti th mi mh e0 da kap Lmin rl1 beta_share 
global T Tmax p q X pop POSN NEGN DN HST Weight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[NP,~] = size(pars);
pars = pars';

%PARAMETERS TO BE ESTIMATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%
bw = pars(1,:).*beta_share;
bs = pars(1,:).*(1-beta_share);
g = pars(2,:);
pp = pars(3,:);
pm = pars(4,:).*pars(3,:);
rl2 = pars(5,:);
c = pars(6,:);

%FIXED PARAMETERS
%%%%%%%%%%%%%%%%%
e0 = e0.*ones(1, NP);
rl1 = rl1.*ones(1, NP);
ti = ti.*ones(1, NP);
th = th.*ones(1, NP);
mi = mi.*ones(1, NP);
mh = mh.*ones(1, NP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = ones(T^2,NP);
alpha2 = ones(T^2,NP);

SA = zeros(T+1,NP);
SR = zeros(T,NP);
S = zeros((T+1)^2,NP);
IA = zeros(T+1,NP);
IR = zeros(T,NP);
I = zeros((T+1)^2,NP);
RA = zeros(T+1,NP);
R = zeros((T+1)^2,NP);
IT = zeros(T+1,NP);
RT = zeros(T+1,NP);
H = zeros(T+1,NP);
tau = zeros(T,NP);

S(1,:) = ones(1,NP).*(1-e0);
I(1,:) = ones(1,NP).*e0;
SA(1,:) = S(1,:);
IA(1,:) = I(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sx = zeros((T+1)*T,NP);
Ix = zeros((T+1)*T,NP);
ITx = zeros(T,NP);
Rx = zeros((T+1)*T,NP);
RTx = zeros(T,NP);
Hx = zeros(T,NP);

Sx((T+1)*(T-1)+1:(T+1)*T,:) = da*(c/2+1)*p(T).*ones(T+1,NP);
Ix((T+1)*(T-1)+1:(T+1)*T-1,:) = da*(c/2+1)*p(T).*ones(T,NP);
ITx(T,:) = da*(c/2+1)*p(T);
Rx((T+1)*(T-1)+1:(T+1)*T-1,:) = da*(c/2+1)*p(T).*ones(T,NP);
RTx(T,:) = da*(c/2+1)*p(T);
Hx(T,:) = da*(c/2+1)*p(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = ones(T,NP);
K0 = 54;
K1 = 70;
K2 = 100;
for t=K0+1:K1
    %L(t,:) = sqrt(rl1.^(t-1-K0)+Lmin^2.*(1-rl1.^(t-1-K0)));
    L(t,:) = sqrt(1 - ((1-Lmin)./(((K1/(K0+1))-1).^rl1)).*(((t/(K0+1))-1).^rl1));
end
for t=K1+1:K2
    L(t,:) = Lmin;
end
for t=K2+1:T
    L(t,:) = sqrt(Lmin^2*rl2.^(t-1-K2)+1-rl2.^(t-1-K2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoOfIterations = 30;
Tolerance = 1e-4;

for k=1:NoOfIterations   
    %Update the State Variables
    for t=1:T
        %Update S variable
        SR(t,:) = sum(alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:),1);
        SA(t+1,:)= SA(t,:)-bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)-bs.*SR(t,:).*IR(t,:);
        %Equation for S^k_{t+1}
        S((T+1)*t+1:(T+1)*t+t,:) = (1-g.*tau(t,:)).*(S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-bw.*(L(t,:)).^2.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IA(t,:)-bs.*alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:));
        %Equation for S^{t+1}_{t+1}
        S((T+1)*t+t+1,:) = g.*tau(t,:).*SA(t+1,:);
        
        %Update I variable
        IR(t,:) = sum(alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:),1);
        IA(t+1,:)= (1-tau(t,:)).*(IA(t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:));
        %Update I^k_{t+1}. Note - I^{t+1}_{t+1} is set to 0.
        I((T+1)*t+1:(T+1)*t+t,:)= (1-tau(t,:)).*(I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IA(t,:)+bs.*alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:));   

        %Update R variable
        RA(t+1,:)= (1-g.*tau(t,:)).*(RA(t,:)+(1-mi)./ti.*IA(t,:));
        %Update R^k_{t+1}. Note - R^{t+1}_{t+1} is set to 0.
        R((T+1)*t+1:(T+1)*t+t,:) = (1-g.*tau(t,:)).*(R((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+(1-mi)./ti.*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:));    

        %Update IT variable
        IT(t+1,:) = (1-1./ti).*IT(t,:)+tau(t,:).*((1-1./ti).*IA(t,:)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:));
        
        %Update RT variable
        RT(t+1,:) = RT(t,:)+(1-mi)./ti.*IT(t,:)+tau(t,:).*g.*(RA(t,:)+(1-mi)./ti.*IA(t,:))+(1-mh)./th.*H(t,:);
        
        %Update H variable
        H(t+1,:) = (1-1./th).*H(t,:)+mi./ti.*(IA(t,:)+IT(t,:));
        
        %Update the testing varible (technically not a state)
        %tau(t,:) = (X(t)-mi./ti.*IA(t,:))./(g.*SA(t+1,:)+g.*(RA(t,:)+(1-mi)./ti.*IA(t,:))+(IA(t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:)));
        tau(t,:) = (X(t))./(g.*SA(t+1,:)+g.*(RA(t,:)+(1-mi)./ti.*IA(t,:))+(IA(t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:)));
    end
    
    %Update the Agent's Adjoint Equations
    for t=1:T-1
        %FOC for da^{t-1}*Sx^k_t
        Sx((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t+1,:) = da*((1-g.*tau(T+1-t,:)).*Sx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)+g.*tau(T+1-t,:).*Sx((T+1)*(T-t)+T+2-t,:))...
            +da*((1-tau(T+1-t,:)).*Ix((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)+tau(T+1-t,:).*ITx(T+1-t,:)-(1-g.*tau(T+1-t,:)).*Sx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)-g.*tau(T+1-t,:).*Sx((T+1)*(T-t)+T+2-t,:)).*(bw.*(L(T+1-t,:)).^2.*IA(T+1-t,:)+bs.*alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:).*IR(T+1-t,:))...
            +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
            -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:)).^2.*(c/2)*q(T+1-t)...
            -da*(1-da)*(L(T+1-t,:)).^2.*IA(T+1-t,:).*bw.*pp*q(T+1-t)...
            -da*(1-da)*alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:).*IR(T+1-t,:).*bs.*pp*q(T+1-t);
        
        %FOC for da^{t-1}*Ix^k_t
        Ix((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t,:) = da*(1-1./ti).*((1-tau(T+1-t,:)).*Ix((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+tau(T+1-t,:).*ITx(T+1-t,:))...
            +da*(1-mi)./ti.*((1-g.*tau(T+1-t,:)).*Rx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+g.*tau(T+1-t,:).*RTx(T+1-t,:))+da*mi./ti.*Hx(T+1-t,:)...
            +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
            -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:)).^2.*(c/2)*q(T+1-t)...
            -da*(1-da)*(L(T+1-t,:)).^2.*SA(T+1-t,:).*bw.*pm*q(T+1-t)...        
            -da*(1-da)*alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:).*SR(T+1-t,:).*bs.*pm*q(T+1-t);
        
        %FOC for da^{t-1}*Rx^k_t
        Rx((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t,:) = da*((1-g.*tau(T+1-t,:)).*Rx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+g.*tau(T+1-t,:).*RTx(T+1-t,:))...
            +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
            -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:)).^2.*(c/2)*q(T+1-t);
        
        %FOC for da^{t-1}*IT_t
        ITx(T-t,:) = da*((1-1./ti).*ITx(T+1-t,:)+(1-mi)./ti.*RTx(T+1-t,:)+mi./ti.*Hx(T+1-t,:))+da*p(T-t)+da*(c/2)*p(T-t);
        
        %FOC for da^{t-1}*RT_t
        RTx(T-t,:) = da*RTx(T+1-t,:)+da*(1-da)*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t);
        
        %FOC for da^{t-1}*H_t
        Hx(T-t,:) = da*((1-1./th).*Hx(T+1-t,:)+(1-mh)./th.*RTx(T+1-t,:))+da*p(T-t)+da*(c/2)*p(T-t);    
    end
 
    %Update the Agent's Control Variable - alpha^k_t
    for t=1:T    
        a = bs./c.*pp.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:)+bs./c.*pm.*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*SR(t,:)+bs./(c.*q(t)*(1-da)).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:).*((1-tau(t,:).*g).*Sx((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+tau(t,:).*g.*Sx((T+1)*(t-1)+t+1,:)-(1-tau(t,:)).*Ix((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-tau(t,:).*ITx(t,:));
        a = 1-a./(1e-12+S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+R((T+1)*(t-1)+1:(T+1)*(t-1)+t,:));  
        alpha2((T)*(t-1)+1:(T)*(t-1)+t,:) = max(0,min(1,a));
    end
    b  = 0.99;
    alpha = b*alpha + (1-b)*alpha2;
    distance = sum(abs(alpha2-alpha),1)/T^2;
    if max(distance) < Tolerance
        break
    end

end

ap2 = 1/(max(POSN)/pop*1e-6)^2;
an2 = 0/(max(NEGN)/pop*1e-6)^2;
ad2 = 1/(max(DN)/pop*1e-6)^2;
ah2 = 0/(max(HST)/pop*1e-6)^2;

XN0 = (RA(2:T+1,:)+SA(2:T+1,:)).*tau.*g;
XP0 = tau(1:T,:).*((1-1./ti).*IA(1:T,:)+bw.*(L(1:T,:)).^2.*SA(1:T,:).*IA(1:T,:)+bs.*SR(1:T,:).*IR(1:T,:))+mi./ti.*IA(1:T,:);
D0 = mh./th.*H(1:T,:);

SS = zeros(1,NP);
% for i=1:NP
% Moment1 = sum(XP0(1:Tmax,:).*(pop*1e6)-POSN)./Tmax;
% Moment2 = sum(XN0(1:Tmax,:).*(pop*1e6)-NEGN)./Tmax;
% Moment3 = sum(D0(1:Tmax,:).*(pop*1e6)-DN)./Tmax;
% Moment4 = sum(H(1:Tmax,:).*(pop*1e6)-HST(1:Tmax))./Tmax;
% Moments = [Moment1; Moment2; Moment3; Moment4];
% SS = sqrt(sum(Moments.*(Weight*Moments)));
% end
SS = SS+sum(((XP0(1:Tmax,:))-((POSN)/(pop*1e6))).^2)*ap2 + sum(((XN0(1:Tmax,:))-((NEGN)/(pop*1e6))).^2)*an2 + sum(((D0(1:Tmax,:))-((DN)/(pop*1e6))).^2)*ad2 + sum(((H(1:Tmax,:))-((HST(1:Tmax))/(pop*1e6))).^2)*ah2;
% SS = SS+sum(((D0(1:Tmax,:))-((DN)/(pop*1e6))).^2)*ad2;
% SS = SS+sum(((H(1:Tmax,:))-((HST(1:Tmax))/(pop*1e6))).^2)*ah2;
SS = SS';

% F = zeros(T+1,T);
% for t=1:T
% F(1:t,t) = ((1-tau(t,:).*g).*Sx((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+tau(t,:).*g.*Sx((T+1)*(t-1)+t+1,:)-(1-tau(t,:)).*Ix((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-tau(t,:).*ITx(t,:));
% end

%FIXED PARAMETERS REDEFINED
%%%%%%%%%%%%%%%%%%%%%%%%%%%
e0 = e0(1, 1);
rl1 = rl1(1, 1);
ti = ti(1, 1);
th = th(1, 1);
mi = mi(1, 1);
mh = mh(1, 1);

end

