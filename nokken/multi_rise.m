clc
clear all
close all

%% Data inladen
out = load('values.mat');
out_zonder= load('values-met-veer-Fv0_0.mat');
T = 0.5;
W = out.w;
m = out.mass;
heffing = out.S.*0.001;
multirise_heffing = [heffing,heffing,heffing,heffing,heffing,heffing,heffing,heffing];
tijd = out.theta/W;
pas  = tijd(2) - tijd(1);
multirise_tijd = [tijd,tijd+0.5,tijd+2*(0.5),tijd+3*(0.5),tijd+4*(0.5),tijd+5*(0.5),tijd+6*(0.5),tijd+7*(0.5)];

% verschil = zeros(length(multirise_tijd),1);
% verschiltijd = zeros(length(multirise_tijd),1);
% xwaarde = zeros(length(multirise_tijd),1);
% for i = 1:(length(multirise_tijd)-2)
%     verschil(i) = multirise_tijd(i+1) - multirise_tijd(i);
%     if i < 36000
%         verschiltijd(i) = tijd(i+1) - tijd(i);
%     end
%     xwaarde(i) = i;
% end
% figure
% hold on
% plot(xwaarde, verschil)
% plot(xwaarde, verschiltijd)
% hold off

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
gamma_numeriek = lsim(sys,theta,tau).';
figure
plot(tau, gamma_numeriek)
ylabel('\gamma_{numeriek} [-]')
xlabel('\tau [-]')

figure
plot(tau, gamma_numeriek- theta)
xlabel('\tau [-]')
ylabel('\gamma_{numeriek}-\theta [-]')

%% Analytische oplossing


% N = 1000;
% x = theta;
% 
% X=fft(x,N);
% a=2*real(X)/N;
% b=-2*imag(X)/N;
% a0 = a(1)/2;
% ak = a(2:N)';
% bk = b(2:N)';
% theta_analytisch = zeros(288000,1);
% for i = 1:288000
%     som = 0;
%     for j = 1:N-1
%         som = som + ak(j)*cos(2*pi*j*tau(i)) + bk(j)*sin(2*pi*j*tau(i));
%     end
%     theta_analytisch(i) = a0 + som;
% end
% figure 
% hold on
% plot(tau, theta_analytisch)
% xlabel('tau')
% %ylabel('theta_{analytisch}')
% % plot(tau, theta)
% % legend('theta_{analytisch}','theta_{numeriek}')
% hold off


K = 100;
% Fourierreeks opstellen,
A = zeros(36000,2*K+1); %288000= length( theta), 201 = 2*K+1 met k= 100 slides 2.27 en 9.39
for i= 1:36000
    for j = 1:2*K+1
        if mod( j , 2 ) == 0
            A(i,j) = cos(2*pi*j*tau(i));
        end
        if mod( j , 2 ) == 1
             if j == 1 
                 A(i,j) = 1/2;
             else
                A(i,j) = sin(2*pi*j*tau(i));
             end
        end
    end
end
X = A\transpose(theta(1:36000));

a_0 = X(1);
coef_a = zeros(K,1);
coef_b = zeros(K,1);
coef_c = zeros(K,1);
coef_d = zeros(K,1);

for i = 2:2*K+1
    if mod(i,2) ==0
        coef_a(i/2) = X(i);
    else
        coef_b((i-1)/2) = X(i);
    end
    
end
for i = 1:K
    coef_c(i) = (-2*zeta*i*coef_b(i)/lambda_tilde+(1-(i/lambda_tilde)^2)*coef_a(i))/((2*zeta*i/lambda_tilde)^2 + (1-(i/lambda_tilde)^2)^2);
    coef_d(i) = (2*zeta*i*coef_a(i)/lambda_tilde+(1-(i/lambda_tilde)^2)*coef_b(i))/((2*zeta*i/lambda_tilde)^2 + (1-(i/lambda_tilde)^2)^2);   
end

gamma_analytisch = zeros(36000,1);
for i = 1:36000
    som = 0;
    for j = 1:K
        som = som + coef_c(j)*cos(2*pi*j*tau(i)) + coef_d(j)*sin(2*pi*j*tau(i));
    end
    gamma_analytisch(i) = a_0/2 + som;
end
theta_analytisch = zeros(36000,1);
for i = 1:36000
    som = 0;
    for j = 1:K
        som = som + coef_a(j)*cos(2*pi*j*tau(i)) + coef_b(j)*sin(2*pi*j*tau(i));
    end
    theta_analytisch(i) = a_0/2 + som;
end
figure 
hold on
plot(tau(1:36000), theta_analytisch)
xlabel('\tau [-]')
%ylabel('theta_{analytisch}')
plot(tau(1:36000), theta(1:36000))
legend('\theta_{analytisch}','\theta_{numeriek}')
hold off


figure
plot(tau(1:36000), gamma_analytisch)
xlabel('\tau [-]')
ylabel('\gamma_{analytisch} [-]')
figure
plot(tau(1:36000), gamma_numeriek(1:36000) - transpose(gamma_analytisch))
xlabel('\tau [-]')
ylabel('\gamma_{numeriek}-\gamma_{analytisch} [-]')


%% Vergelijken multi- en singlerise
sr = load('single_rise.mat');
gamma_sr = sr.gamma_numeriek;

% We kiezen periode 4
gamma_mr = gamma_numeriek(128000:136000);
theta_mr = theta(128000:136000);
tau4 = tau(128000:136000);

figure
hold on
plot(tau4, gamma_mr, 'b')
plot(tau4, 1-gamma_sr, 'r')
ylabel('\gamma [-]')
xlabel('\tau [-]')
legend('\gamma_{mr}', '\gamma_{sr}');
hold off

figure 
plot(tau4, (gamma_mr-(1-gamma_sr)))
ylabel('\gamma_{mr}-\gamma_{sr} [-]')
xlabel('\tau [-]')

%meer ingezoomde figuren om enkel vrije respons te vergelijken
% single rise -> 8000 waarden, tau = 1 voor index 4001
% multi rise periode 4 -> waarden 132001 tot 136000
gamma_sr_vrijerespons = sr.gamma_numeriek(4001:8000);
gamma_mr_vrijerespons = gamma_numeriek(132001:136000);
theta_mr_vrijerespons = theta(132001:136000);
tau4_vrijerespons = tau(132001:136000);

figure
hold on
plot(tau4_vrijerespons, gamma_mr_vrijerespons, 'b')
ylabel('\gamma_{mr-vrije respons} [-]')
xlabel('\tau [-]')

hold off

figure 
plot(tau4_vrijerespons, (gamma_mr_vrijerespons-(1-gamma_sr_vrijerespons)))
ylabel('\gamma_{mr-vrije respons}-\gamma_{sr-vrije respons} [-]')
xlabel('\tau [-]')
%% nieuwe veerkracht met vervormbare volger
F = zeros(length(out_zonder.theta),1);
gamma_veer= gamma_numeriek(1:36000);
theta_veer= theta(1:36000);
for i= 1:length(out_zonder.theta)
    F(i) = (0*out_zonder.S(i) - (gamma_veer(i) - theta_veer(i))*kf*0.001*40)/cos(out_zonder.pressure_angle(i));
    Fm(i)= (200 + 275*out_zonder.S(i) - (gamma_veer(i) - theta_veer(i))*kf*0.001*40)/cos(out_zonder.pressure_angle(i));
end


figure
hold on
plot(out_zonder.thetadegree,F)
plot(out_zonder.thetadegree,Fm)
hold off
%plot(out_zonder.thetadegree, out_zonder.normalforce_tot)
legend('F_{zonder veer}', 'F_{met veer}');
xlabel('theta [-]')
ylabel('F van mij')

