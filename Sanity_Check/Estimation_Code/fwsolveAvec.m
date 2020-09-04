function [SS,XP0,XN0,D0,SA,SR,L, tau, alpha] = fwsolveAvec(pars)
global da T Tmax p q X Lmin pop POSN DN HST

[NP,~] = size(pars);
pars = pars';

eta = pars(1,:);
bw = pars(2,:).*eta;
bs = pars(2,:).*(1-eta);
g = pars(3,:);
pp = pars(4,:);
pm = pars(5,:).*pars(4,:);
rl1 = pars(6,:);
rl2 = pars(7,:);
e0 = pars(8,:)/pop*1e-6;
c = pars(9,:);
ti = pars(10,:);
th = pars(11,:);
mi = pars(12,:);
mh = pars(13,:);

alpha = ones(T^2,NP);
alpha2 = ones(T^2,NP);

S = zeros((T+1)^2,NP);
I = zeros((T+1)^2,NP);
IT = zeros(T+1,NP);
H = zeros(T+1,NP);
SA = zeros(T+1,NP);
IA = zeros(T+1,NP);
SR = zeros(T,NP);
IR = zeros(T,NP);
R = zeros((T+1)^2,NP);
RT = zeros(T+1,NP);
RA = zeros(T+1,NP);
tau = zeros(T,NP);

Sx = zeros((T+1)*T,NP);
Sx((T+1)*(T-1)+1:(T+1)*T,:) = da*(c/2+1)*p(T).*ones(T+1,NP);
Ix = zeros((T+1)*T,NP);
Ix((T+1)*(T-1)+1:(T+1)*T-1,:) = da*(c/2+1)*p(T).*ones(T,NP);
ITx = zeros(T,NP);
ITx(T,:) = da*(c/2+1)*p(T);
Rx = zeros((T+1)*T,NP);
Rx((T+1)*(T-1)+1:(T+1)*T-1,:) = da*(c/2+1)*p(T).*ones(T,NP);
RTx = zeros(T,NP);
RTx(T,:) = da*(c/2+1)*p(T);
Hx = zeros(T,NP);
Hx(T,:) = da*(c/2+1)*p(T);

S(1,:) = ones(1,NP).*(1-e0);
I(1,:) = ones(1,NP).*e0;
SA(1,:) = S(1,:);
IA(1,:) = I(1,:);

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



for k=1:30
    
for t=1:T
    SR(t,:) = sum(alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:),1);
    IR(t,:) = sum(alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:),1);
    SA(t+1,:)= SA(t,:)-bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)-bs.*SR(t,:).*IR(t,:);
    tau(t,:) = (X(t)-mi./ti.*IA(t,:))./(g.*SA(t+1,:)+g.*(RA(t,:)+(1-mi)./ti.*IA(t,:))+(IA(t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:)));
    RA(t+1,:)= (1-g.*tau(t,:)).*(RA(t,:)+(1-mi)./ti.*IA(t,:));
    IA(t+1,:)= (1-tau(t,:)).*(IA(t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:));
    S((T+1)*t+t+1,:) = g.*tau(t,:).*SA(t+1,:);
    IT(t+1,:) = (1-1./ti).*IT(t,:)+tau(t,:).*((1-1./ti).*IA(t,:)+bw.*(L(t,:)).^2.*SA(t,:).*IA(t,:)+bs.*SR(t,:).*IR(t,:));     
    RT(t+1,:) = RT(t,:)+(1-mi)./ti.*IT(t,:)+tau(t,:).*g.*(RA(t,:)+(1-mi)./ti.*IA(t,:))+(1-mh)./th.*H(t,:);
    H(t+1,:) = (1-1./th).*H(t,:)+mi./ti.*(IA(t,:)+IT(t,:));
    S((T+1)*t+1:(T+1)*t+t,:) = (1-g.*tau(t,:)).*(S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-bw.*(L(t,:)).^2.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IA(t,:)-bs.*alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:));
    I((T+1)*t+1:(T+1)*t+t,:)= (1-tau(t,:)).*(I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*(1-1./ti)+bw.*(L(t,:)).^2.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IA(t,:)+bs.*alpha((T)*(t-1)+1:(T)*(t-1)+t,:).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:));   
    R((T+1)*t+1:(T+1)*t+t,:) = (1-g.*tau(t,:)).*(R((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+(1-mi)./ti.*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:));    
end

for t=1:T-1
    Sx((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t+1,:) = da*((1-g.*tau(T+1-t,:)).*Sx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)+g.*tau(T+1-t,:).*Sx((T+1)*(T-t)+T+2-t,:))...
        +da*((1-tau(T+1-t,:)).*Ix((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)+tau(T+1-t,:).*ITx(T+1-t,:)-(1-g.*tau(T+1-t,:)).*Sx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t+1,:)-g.*tau(T+1-t,:).*Sx((T+1)*(T-t)+T+2-t,:)).*(bw.*(L(T+1-t,:)).^2.*IA(T+1-t,:)+bs.*alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:).*IR(T+1-t,:))...
        +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:)).^2.*(c/2)*q(T+1-t)...
        -da*(1-da)*(L(T+1-t,:)).^2.*IA(T+1-t,:).*bw.*pp*q(T+1-t)...
        -da*(1-da)*alpha((T)*(T-t)+1:(T)*(T-t)+T-t+1,:).*IR(T+1-t,:).*bs.*pp*q(T+1-t);
    Ix((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t,:) = da*(1-1./ti).*((1-tau(T+1-t,:)).*Ix((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+tau(T+1-t,:).*ITx(T+1-t,:))...
        +da*(1-mi)./ti.*((1-g.*tau(T+1-t,:)).*Rx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+g.*tau(T+1-t,:).*RTx(T+1-t,:))+da*mi./ti.*Hx(T+1-t,:)...
        +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:)).^2.*(c/2)*q(T+1-t)...
        -da*(1-da)*(L(T+1-t,:)).^2.*SA(T+1-t,:).*bw.*pm*q(T+1-t)...        
        -da*(1-da)*alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:).*SR(T+1-t,:).*bs.*pm*q(T+1-t);
    Rx((T+1)*(T-t-1)+1:(T+1)*(T-t-1)+T-t,:) = da*((1-g.*tau(T+1-t,:)).*Rx((T+1)*(T-t)+1:(T+1)*(T-t)+T-t,:)+g.*tau(T+1-t,:).*RTx(T+1-t,:))...
        +da*(1-da).*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t)...
        -da*(1-da)*(1-alpha((T)*(T-t)+1:(T)*(T-t)+T-t,:)).^2.*(c/2)*q(T+1-t);
    ITx(T-t,:) = da*((1-1./ti).*ITx(T+1-t,:)+(1-mi)./ti.*RTx(T+1-t,:)+mi./ti.*Hx(T+1-t,:))+da*p(T-t)+da*(c/2)*p(T-t);
    RTx(T-t,:) = da*RTx(T+1-t,:)+da*(1-da)*L(T+1-t,:)*q(T+1-t)+da*p(T-t)+da*(1-da)*(c/2)*q(T+1-t)+da*(c/2)*p(T-t);
    Hx(T-t,:) = da*((1-1./th).*Hx(T+1-t,:)+(1-mh)./th.*RTx(T+1-t,:))+da*p(T-t)+da*(c/2)*p(T-t);
    
end
 
for t=1:T    
        a = bs./c.*pp.*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:)+bs./c.*pm.*I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*SR(t,:)+bs./c./(q(t)*(1-da)).*S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:).*IR(t,:).*((1-tau(t,:).*g).*Sx((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+tau(t,:).*g.*Sx((T+1)*(t-1)+t+1,:)-(1-tau(t,:)).*Ix((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-tau(t,:).*ITx(t,:));
        a = 1-a./(1e-24+S((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+I((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+R((T+1)*(t-1)+1:(T+1)*(t-1)+t,:));  
        alpha2((T)*(t-1)+1:(T)*(t-1)+t,:) = max(0,min(1,a));
end

diffa = sum(abs(alpha2-alpha),1)/T^2;
b  = 0.9;
alpha = b*alpha + (1-b)*alpha2;

if max(diffa) < 1e-4
    break
end

end

ap2 = 1/(max(POSN(1:Tmax))/pop*1e-6)^2;
ad2 = 1/(max(DN(1:Tmax))/pop*1e-6)^2;
ah2 = 0/(max(HST(1:Tmax))/pop*1e-6)^2;

XP0 = tau(1:T,:).*((1-1./ti).*IA(1:T,:)+bw.*(L(1:T,:)).^2.*SA(1:T,:).*IA(1:T,:)+bs.*SR(1:T,:).*IR(1:T,:))+mi./ti.*IA(1:T,:);
XN0 = (RA(2:T+1,:)+SA(2:T+1,:)).*tau.*g;
D0 = mh./th.*H(1:T,:);

SS = zeros(1,NP);


SS = SS+sum(((XP0(1:Tmax,:))-(POSN(1:Tmax))/pop*1e-6).^2)*ap2;
SS = SS+sum(((D0(1:Tmax,:))-(DN(1:Tmax))/pop*1e-6).^2)*ad2;
SS = SS+sum((H(1:Tmax,:)-HST(1:Tmax)/pop*1e-6).^2)*ah2;

SS = SS';

% F = zeros(T+1,T);
% for t=1:T
% F(1:t,t) = ((1-tau(t,:).*g).*Sx((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)+tau(t,:).*g.*Sx((T+1)*(t-1)+t+1,:)-(1-tau(t,:)).*Ix((T+1)*(t-1)+1:(T+1)*(t-1)+t,:)-tau(t,:).*ITx(t,:));
% end

end

