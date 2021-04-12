out = load('values.mat');
out_exc = load('values_exc.mat');
P_exc = zeros(length(out.theta),1);
P = zeros(length(out.theta),1);
P_average = zeros(length(out.theta),1);
P_average_exc = zeros(length(out.theta),1);
M_last = zeros(length(out.theta),1);
P_tot_exc = 0;
P_tot = 0;
for i = 1:length(out.theta)
   
    P(i) = (10^-3)*out.normalforce_tot(i)*sin(out.pressure_angle(i))*(out.bcr+ out.rof+out.S(i))*out.w;
    P_tot = P_tot + P(i);
    M_last(i) = P(i)/out.w;
    P_exc(i) = (10^-3)*out_exc.normalforce_tot(i) * (out_exc.exc*cos(out_exc.pressure_angle(i))+sin(out_exc.pressure_angle(i))*(sqrt((out_exc.bcr + out_exc.rof)^2-(out_exc.exc)^2)+out_exc.S(i)))*out_exc.w;
    P_tot_exc = P_tot_exc + P_exc(i);
end


P_average_exc = P_average_exc + P_tot_exc/length(out.theta);
P_average = P_average + P_tot/length(out.theta);

% 
% figure
% plot(theta, P_exc)
% xlabel('theta')
% ylabel('P_{exc}')
% figure
% plot(theta, P)
% xlabel('theta')
% ylabel('P')
% 
% figure
% plot(theta, P-P_exc)
% xlabel('theta')
% ylabel('P-P_{exc}')
% 
figure
plot(out.theta, M_last, 'g')
hold on
plot(out.theta,P_average/out.w, 'r')
hold off
xlabel('crank angle [°]')

legend('M_{last}','M_{drive}')


A = zeros(length(out.theta),2);
A = cumtrapz(out.theta,M_last-P_average/out.w);

figure
plot(out.theta, A)
ylabel('A[J]')
xlabel('crank angle [°]')

theta_max = 3.075;
theta_min = 1.777;
index_max = find(abs(out.theta - theta_max)< 2*pi/36000 );
index_min = find(abs(out.theta - theta_min)< 2*pi/36000 );

A_max = trapz(out.theta(index_min:index_max), M_last(index_min:index_max)-P_average(index_min:index_max)/out.w);
K = 0.1;
I = A_max/(K*out.w^2); %w_av = out.w, eenheid I = [kg*m^2]

%kies R = 225mm = basiscirkel + roller radius + max lift
R = 0.225 ;
rho = 7800; 
t = 2*I/(pi*R^4*rho);


