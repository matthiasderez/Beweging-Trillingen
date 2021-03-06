clc
clear all
close all

%% Data inladen
out = load('values.mat');

%selecteren van de waardes bij de fall tss hoek 200? en 240?
heffing = out.S(20000:28000)*0.001;
tijd = out.theta(20000:28000)/out.w;
m = out.mass;


%% Variabelen aanmaken
t1 = 40/720; %zie andreas
T = (2*pi)/out.w;
tau = (tijd - tijd(1))/t1; % -tijd(1) om ervoor te zorgen dat ons segment start op tau = 0
step = tau(2)-tau(1);
zeta = 0.091; %gegeven
lambda = 0.75/zeta; % 10% accuraat
lambda_d = lambda*sqrt(1-zeta^2);
theta = (0.04-heffing)./0.04; %van fall een rise maken

%% kf berekenen
kf = m*(lambda*2*pi/t1)^2;
%% Transferfunctie aanmaken
teller = (2*pi*lambda)^2;
noemer = [1, 2*zeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(teller, noemer);

%% Plotten figuren 
figure
plot(tau, theta)
xlabel('\tau')
ylabel('\theta')


%% Numeriek oplossen met lsim
figure
lsim(sys, theta, tau)
gamma_numeriek = lsim(sys, theta, tau);
gamma_numeriek = transpose(gamma_numeriek);
figure
plot(tau, (gamma_numeriek-theta));
xlabel('\tau [-]');
ylabel('\gamma_{numeriek}-\theta [-]')
figure
plot(tau, gamma_numeriek)
xlabel('\tau [-]')
ylabel('\gamma_{numeriek} [-]')
figure
plot(tau, gamma_numeriek)
xlim([1 2])
xlabel('\tau [-]')
ylabel('\gamma_{numeriek} [-]')



%% Analystisch oplossen 
% voor vrije respons dus tau > 1
% zie slides 12 


[value,index]=min(abs(tau-1)); % index vinden waar tau zo dicht mogelijk ligt bij 1
gamma1 = gamma_numeriek(index);
gammadot1 = (gamma_numeriek(index+1)-gamma_numeriek(index-1))./(2*step); %afgeleide numeriek benaderen
%formules slide13
A = sqrt((((gamma1-1)*2*pi*lambda_d)^2+(gammadot1+zeta*2*pi*lambda*(gamma1-1))^2)/(2*pi*lambda_d)^2);
phi = atan(-(gammadot1+zeta*2*pi*lambda*(gamma1-1))/((gamma1-1)*2*pi*lambda_d));

oscillatie = zeros(1,length(tau));
compl_omh = zeros(1,length(tau));
for i = 1:length(tau)
    oscillatie(i) = A*exp(-zeta*2*pi*lambda*(tau(i)-1))*cos(2*pi*lambda_d*(tau(i)-1) + phi);
    compl_omh(i) = A*exp(-zeta*2*pi*lambda*(tau(i)-1));
end
gamma_analytisch = 1+ oscillatie;
figure
hold on
plot(tau, gamma_analytisch)
xlim([1 2])
xlabel('\tau [-]')
ylabel('\gamma_{analytisch} [-]')
plot(tau, 1+compl_omh)
plot(tau, 1-compl_omh)
hold off

figure
plot(tau(4001:8001), (gamma_analytisch(4001:8001) - gamma_numeriek(4001:8001))./gamma_numeriek(4001:8001))
xlim([1 2]);
xlabel('\tau [-]')
ylabel('(\gamma_{analytisch} - \gamma_{numeriek})/ \gamma_{numeriek})[-]')

%% Benaderende oplossing
%Voorwaarde checken
controle = exp(-zeta*2*pi*lambda);
Q = (2*pi)^2;
N = 3;
Ab = Q/(2*pi*lambda)^N * sqrt(1/(1-zeta^2)); %formule slide 27
epsilon = abs((A-Ab)/A);

b2 = -zeta*2*pi/lambda;
b1 = -1/lambda^2*(1-4*zeta^2);
b0 = 2*zeta/(pi*lambda^3)*(1-2*zeta^2)+1;  


theta_0 = b0-1; %zie notities matthias + slide 9.27 theta_0 = gamma(1) - 1
theta_dot0 = b1; %zie notities matthias + slide 9.27 theta_dot0 = gammadot(1) 
[A_appr,B_appr,C_appr,D_appr] = tf2ss(teller,noemer);
X0 = [1/C_appr(2)*theta_dot0; 1/C_appr(2)*theta_0];
gamma_b = transpose(lsim(A_appr,B_appr,C_appr,D_appr, theta, tau, X0));
figure
plot(tau, gamma_b)
xlim([1 2]);
xlabel('\tau [-]')
ylabel('\gamma_{benaderend} [-]')

figure
plot(tau, (gamma_b - gamma_numeriek)./gamma_numeriek)
xlim([1 2]);
xlabel('\tau [-]')
ylabel('(\gamma_{benaderend}-\gamma_{numeriek})/\gamma_{numeriek} [-]')





%%% Numeriek op andere manier gedaan, fall niet als rise beschouwd%%%


% theta = heffing./40;
% theta0 = 1; % vul hier zelf de initiele dimensieloze heffing in
% theta_dot0 = 0; % vul hier zelf de initiele dimensieloze snelheid in
% [A1,B,C,D] = tf2ss(teller,noemer);
% X0 = [1/C(2)*theta_dot0; 1/C(2)*theta0];
% lsim(A1,B,C,D, theta, tau, X0)
% xlabel('tau (-)')
% ylabel('gamma (tau)')
% gamma = lsim(A1,B,C,D, theta, tau, X0);
% 
% figure
% plot(tau, transpose(gamma)-theta)
% xlabel('tau(-)');
% ylabel('gamma-theta')
% 
% 
% 
% 
% 
% index = find(tau==1);
% gamma1 = gamma(index);
% 
% gammadot1 = (gamma(index+1)-gamma(index-1))./(2*step);
% lambda_d = lambda*sqrt(1-zeta^2);
% A2 = sqrt((((gamma1)*2*pi*lambda_d)^2+(gammadot1+zeta*2*pi*lambda*(gamma1))^2)/(2*pi*lambda_d)^2);
% %opm in de formule erboven (slide 9.13) x_0 = gamma1 ipv gamma1-1 omdat wij
% %zitten oscilleren rond nul ipv 1
% phi = atan(-(gammadot1+zeta*2*pi*lambda*gamma1)/(gamma1*2*pi*lambda_d)) + pi;
% 
% oscillatie2 = zeros(length(tau),1);
% compl_omh2 = zeros(length(tau),1);
% for i = 1:length(tau)
%     oscillatie2(i) = A*exp(-zeta*2*pi*lambda*(tau(i)-1))*cos(2*pi*lambda_d*(tau(i)-1) + phi);
%     compl_omh2(i) = A*exp(-zeta*2*pi*lambda*(tau(i)-1));
% end
% 
% figure
% hold on
% plot(tau, oscillatie)
% plot(tau, compl_omh)
% plot(tau, -compl_omh)
% hold off
% 
% 
save('single_rise.mat')