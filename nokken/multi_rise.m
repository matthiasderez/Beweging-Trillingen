clc
clear all
close all

%% Data inladen
out = load('values.mat');
T = 0.5;
W = out.w;
m = out.mass;
heffing = out.S.*0.001;
multirise_heffing = [heffing,heffing,heffing,heffing,heffing,heffing,heffing,heffing];
tijd = out.theta/W;
tijd
tijd+1
multirise_tijd = [tijd,tijd+0.5/36000 + 0.5,tijd+0.5/36000 + 2*0.5,tijd+0.5/36000 + 3*0.5,tijd+0.5/36000 + 4*0.5,tijd+0.5/36000 +5* 0.5,tijd+0.5/36000 + 6*0.5,tijd+0.5/36000 +7* 0.5];

%% Dedimensionaliseren
theta = multirise_heffing/0.04;
tau = multirise_tijd/T;

%% Variabelen
t1_sr = 40/(W*180/pi);
zeta = 0.091;
lambda = 0.75/zeta;
kf = m*(lambda*2*pi/t1_sr)^2;
w_n = sqrt(kf/m);
t_n = 2*pi/w_n;
lambda_tilde = T/t_n;

%% Overdrachtsfunctie aanmaken
teller = (2*pi*lambda_tilde)^2;
noemer = [1, 2*zeta*(2*pi*lambda_tilde), (2*pi*lambda_tilde)^2];
sys = tf(teller, noemer);

%% Simulatie
gamma = lsim(sys,theta,tau).';
figure
plot(tau, gamma)
ylabel('gamma')
xlabel('tau')

figure
plot(tau, gamma- theta)
xlabel('tau')
ylabel('gamma-tau')


% 
% figure 
% plot(tau,heffing)
% xlabel('tau[-]')
% ylabel('theta(tau) [-]')


