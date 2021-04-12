clc
clear all
close all

%% import data
datamat = load('cam_design.mat');
W = datamat.w;
T = (2*pi)/W;
m = datamat.mass;
lift = datamat.S*0.001;
extendedlift = [lift,lift,lift,lift,lift,lift,lift,lift,lift,lift,lift];
time = datamat.theta/W;
extendedtime = [time,time+time(36000)+0.1/36000,time+2*(time(36000)+0.1/36000),time+3*(time(36000)+0.1/36000),time+(0.1/36000+time(36000))*4,time+(time(36000)+0.1/36000)*5,time+(time(36000)+0.1/36000)*6,time+(time(36000)+0.1/36000)*7,time+(time(36000)+0.1/36000)*8,time+(time(36000)+0.1/36000)*9,time+(time(36000)+0.1/36000)*10];
%% dedimensionaliseren
theta = extendedlift/0.04;
tau = extendedtime/T;
%% overdrachtsfunctie
zeta = 0.096;
lambdasr = 0.75/zeta;
t1sr = T*(50/360);
kf = ((2*pi*sqrt(m)*lambdasr)/(t1sr))^2;
wn = sqrt(kf/m);
tn = (2*pi)/wn;
lambda = T/tn;
teller = (2*pi*lambda)^2;
noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(teller, noemer);
%% simulatie
gamma = lsim(sys,theta,tau).';
figure
plot(tau, gamma)
ylabel('gamma')
xlabel('tau')
%% neem de 11 periode
gamma11 = gamma(360001:396000);
theta11 = theta(360001:396000);
tau11 = tau(360001:396000);
%% plot 1 periode
figure
plot(tau11, gamma11, 'b')
ylabel('gamma')
xlabel('tau')
figure
plot(tau11, gamma11, 'b',tau11, theta11, 'r')
ylabel('gamma&theta')
xlabel('tau')
figure
plot(tau11, gamma11-theta11)
ylabel('gamma-theta')
xlabel('tau')
%% single rise
liftsr = datamat.S(10000:16000)*0.001;
thetasr = (liftsr-0.02)/0.02;
timesr = datamat.theta(10000:16000)/datamat.w;
tsr =  timesr-timesr(1);
t1sr = T*(50/360);
tausr = tsr/t1sr;
tellersr = (2*pi*lambdasr)^2;
noemersr = [1, 2*zeta*(2*pi*lambdasr), (2*pi*lambdasr)^2];
theta0 = 1; % vul hier zelf de initiele dimensieloze heffing in
theta_dot0 = 0; % vul hier zelf de initiele dimensieloze snelheid in
[A,B,C,D] = tf2ss(tellersr,noemersr);
X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];
gammasr = (lsim(A,B,C,D, thetasr, tausr, X0).'+1)/2;
%% deel dat overeenkomt met sr
gammamr = gamma11(10000:16000);
thetamr = theta11(10000:16000);
%% vergelijking
figure
plot(tausr, gammamr, 'b', tausr, gammasr, 'r')
ylabel('gamma')
xlabel('tau')
figure
plot(tausr, gammamr-thetamr, 'b', tausr, gammasr-thetamr, 'r')
ylabel('gamma-theta')
xlabel('tau')
figure
plot(tausr, gammamr-gammasr)
ylabel('multirise - single rise')
xlabel('tau')
%% Invloed op contact kracht 
% beweging is gamma11
% serieschakeling van kf en k gevonden in deel 3
%k3 = datamat.springconstant;
%kfmm = kf/1000;
%keq = (1/k3+1/kfmm)^(-1);
%verplaatsing = gamma11*40;
%calpha = cos(datamat.pressure_angle);
%Fve = verplaatsing*keq;
%Fvo = datamat.springpreload;
%Fveer = (Fve+Fvo)./calpha;
%Facc = datamat.normalforce_acc;
%Ff = datamat.normalforce_load;
%Ftot = Fveer+Facc+Ff;
%figure
%plot(tau11, Ff ,'g',tau11, Facc ,'r',tau11, Fveer ,'y', tau11, Ftot, 'b')
%ylabel('Normal forces (N)')
%xlabel('tau')
%% Invloed op kracht
% We hebben de functionele krachten, de veerkrachten berekend in deel 3
% en nu ook de kracht door de veer kf
calpha = cos(datamat.pressure_angle);
k3 = datamat.springconstant;
F3 = datamat.springpreload;
Fnv3 = (F3+k3*theta11*40)./calpha;
Fnload = datamat.normalforce_load;
indrukveerf = (theta11-gamma11)*0.040;
Fvkf = (indrukveerf*kf)./calpha;
Ftot = Fvkf + Fnload + Fnv3;
figure
plot(tau11, Fnv3 ,'y',tau11, Fnload ,'r',tau11, Fvkf,'g',tau11, Ftot,'b')
ylabel('Normal forces (N)')
xlabel('tau')
%% Aanpassing k3 en F3
k3nieuw = k3+72;
F3nieuw = F3+900;
Fnv3nieuw = (F3nieuw+k3nieuw*theta11*40)./calpha;
Ftotnieuw = Fvkf + Fnload + Fnv3nieuw;
figure
plot(tau11, Fnv3nieuw ,'y',tau11, Fnload ,'r',tau11, Fvkf,'g',tau11, Ftotnieuw,'b')
ylabel('Normal forces (N)')
xlabel('tau')







