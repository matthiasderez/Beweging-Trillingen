clc
clear all
close all

%% parameters
Q = (2*pi)^2;
N = 3;
zeta = 0.096;
%tbegin = 0;
%teind = 0.1*(60/360);
%tstap = 0.0001;
%t=[tbegin:tstap:teind];

%% benadering minstens 10% accuraat en zo lambda vinden

lambda = 0.75/zeta;
%% importen matcam
datamat = load('cam_design.mat');
lift = datamat.S(10000:16000)*0.001;
time = datamat.theta(10000:16000)/datamat.w;
T = (2*pi)/datamat.w;
m = datamat.mass;
t =  time-time(1);

%% dimensieloze parameters definieren
t1 = T*(50/360);
tau = t/t1;
%% kf uitrekenen
kf = ((2*pi*sqrt(m)*lambda)/(t1))^2/1000;


%% overdrachtsfucntie aanmaken

teller = (2*pi*lambda)^2;
noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(teller, noemer);

%% theta
theta = (lift-0.02)/0.02;
%% simulatie met beginvoorwaarden
theta0 = 1; % vul hier zelf de initiele dimensieloze heffing in
theta_dot0 = 0; % vul hier zelf de initiele dimensieloze snelheid in
[A,B,C,D] = tf2ss(teller,noemer);
X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];
gamma = lsim(A,B,C,D, theta, tau, X0);
gammareken = gamma.';
plotter = gammareken-theta;
figure
plot(tau, gammareken)
ylabel('gamma')
xlabel('tau')
figure
plot(tau, plotter)
ylabel('gamma-theta')
xlabel('tau')
%% vergelijken met benadering

g1 = gamma(5001);
lambdad = lambda*sqrt(1-zeta^2);
x0 = g1;
v0 = (gamma(5002)-gamma(5000))/(tau(5002)-tau(5000));
A1 = sqrt(((x0*2*pi*lambdad)^2+(v0+zeta*2*pi*lambda*x0)^2)/((2*pi*lambdad)^2));

Atilde = (Q/(2*pi*lambda)^3)*sqrt(1/(1-zeta^2));
epsilon = abs((A1-Atilde)/A1);
%% figuren make
gammana1 = gamma(5001:6001);
tauna1 = tau(5001:6001);
omhullende1 = A1*exp(-zeta*2*pi*lambda*(tauna1-1));
omhullendebenadering1 = Atilde*exp(-zeta*2*pi*lambda*(tauna1-1));
omhullende2 = -A1*exp(-zeta*2*pi*lambda*(tauna1-1));
omhullendebenadering2 = -Atilde*exp(-zeta*2*pi*lambda*(tauna1-1));
figure
plot(tauna1, gammana1, 'b', tauna1, omhullende1, 'b', tauna1, omhullende2, 'b')
ylabel('gamma(tau)')
xlabel('tau')

figure
plot(tauna1, gammana1)
ylabel('gamma(tau)')
xlabel('tau')

%% simulatie benadering
%thetab = -(Q/6)*(tau-1).^3;
%thetab0 = 1; % vul hier zelf de initiele dimensieloze heffing in
%theta_dotb0 = 0; % vul hier zelf de initiele dimensieloze snelheid in
%[Ab,Bb,Cb,Db] = tf2ss(teller,noemer);
%Xb0 = [1/Cb(2)*theta_dotb0; 1/Cb(2)*thetab0];
%lsim(Ab,Bb,Cb,Db, thetab, tau, Xb0);
%% freeresponsebenadering1
freeresponse = zeros(1001,1);
taufree = tauna1;
thetabb0 = -(2*Q*zeta)/((2*pi*lambda)^3)*(1-(4*zeta^2-1)); % vul hier zelf de initiele dimensieloze heffing in
theta_dotbb0 = -Q*(4*zeta^2-1)/((2*pi*lambda)^2); % vul hier zelf de initiele dimensieloze snelheid in
[Abb,Bbb,Cbb,Dbb] = tf2ss(teller,noemer);
Xbb0 = [1/Cbb(2)*theta_dotbb0; 1/Cbb(2)*thetabb0];
gammafree = lsim(Abb,Bbb,Cbb,Dbb, freeresponse, taufree, Xbb0);
figure
plot(tauna1, gammafree, 'b', tauna1, omhullendebenadering1, 'b', tauna1, omhullendebenadering2, 'b')
ylabel('gamma approximation(tau)')
xlabel('tau')
figure
plot(tauna1, gammafree, 'b')
ylabel('gamma approximation(tau)')
xlabel('tau')
%% freeresponsebenadering2
%freeresponse = zeros(1001,1);
%taufree = tauna1-1;
%thetabb0 = gammab(5001)-1 % vul hier zelf de initiele dimensieloze heffing in
%theta_dotbb0 = (gammab(5002)-gammab(5000))/(tau(5002)-tau(5000)) % vul hier zelf de initiele dimensieloze snelheid in
%[Abb,Bbb,Cbb,Dbb] = tf2ss(teller,noemer);
%Xbb0 = [1/Cbb(2)*theta_dotbb0; 1/Cbb(2)*thetabb0];
%lsim(Abb,Bbb,Cbb,Dbb, freeresponse, taufree, Xbb0)
%% verschil
diff= gammana1-gammafree;
figure
plot(tauna1,gammana1,'b', tauna1, omhullende1, 'b', tauna1, omhullende2, 'b',tauna1,gammafree,'r', tauna1, omhullendebenadering1, 'r', tauna1, omhullendebenadering2, 'r')
ylabel('free response')
xlabel('tau')
figure
plot(tauna1,diff)
ylabel('numerical-approximate')
xlabel('tau')













