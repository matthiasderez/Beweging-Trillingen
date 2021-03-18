%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi3,phi4,phi5,phi6,phi7,phi8,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8] = kinematics_4bar(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,phi1,phi2,dphi2,ddphi2,phi3_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,t,fig_kin_4bar)


% allocation of the result vectors (this results in better performance because we don't have to reallocate and
% copy the vector each time we add an element.
phi3 = zeros(size(t));
phi4 = zeros(size(t));
phi5 = zeros(size(t));
phi6 = zeros(size(t));
phi7 = zeros(size(t));
phi8 = zeros(size(t));
dphi3 = zeros(size(t));
dphi4 = zeros(size(t));
dphi5 = zeros(size(t));
dphi6 = zeros(size(t));
dphi7 = zeros(size(t));
dphi8 = zeros(size(t));
ddphi3 = zeros(size(t));
ddphi4 = zeros(size(t));
ddphi5 = zeros(size(t));
ddphi6 = zeros(size(t));
ddphi7 = zeros(size(t));
ddphi8 = zeros(size(t));

% fsolve options (help fsolve, help optimset)
optim_options = optimset('Display','off', 'TolFun', 1e-12);

% *** loop over positions ***
Ts = t(2) - t(1);      % timestep
t_size = size(t,1);    % number of simulation steps
for k=1:t_size
    
    % *** position analysis ***
    
    % fsolve solves the non-linear set of equations
    % loop closure equations: see loop_closure_eqs.m
    % argument loop_closure_eqs: file containing closure equations
    % argument [..]': initial values of unknown angles phi3 and phi4
    % argument optim options: parameters for fsolve
    % argument phi2(k): input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
    % argument a1 ... phi1: constants
    % return value x: solution for the unknown angles phi3 and phi4
    % return exitflag: indicates convergence of algorithm
    [x, fval, exitflag]=fsolve('loop_closure_eqs',[phi3_init phi4_init phi5_init phi6_init phi7_init phi8_init]',optim_options,phi2(k),r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,phi1);
    if (exitflag ~= 1)
        display 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    %else
        %display 'convergence!'
    end
    
    % save results of fsolve
    phi3(k)=x(1);
    phi4(k)=x(2);
    phi5(k)=x(3);
    phi6(k)=x(4);
    phi7(k)=x(5);
    phi8(k)=x(6);
    
    
    
    % *** velocity analysis ***
    
    Bm = [-r3*sin(phi3(k)),  r4*sin(phi4(k)),0,0,0,0;
         r3*cos(phi3(k)), -r4*cos(phi4(k)),0,0,0,0;
         0,0,-r5*sin(phi5(k)), -r6*sin(phi6(k)), r7*sin(phi7(k)),-r8*sin(phi8(k));
         0,0,r5*cos(phi5(k)), r6*cos(phi6(k)), -r7*cos(phi7(k)),r8*cos(phi8(k));
         -r3*sin(phi3(k)),0,0, r6*sin(phi6(k)), r10*sin(phi7(k)),0;
         r3*cos(phi3(k)),0,0, -r6*cos(phi6(k)), -r10*cos(phi7(k)),0];
    Am = [ r2*sin(phi2(k))*dphi2(k);
         -r2*cos(phi2(k))*dphi2(k);
         0;
         0;
         -r9*sin(phi2(k))*dphi2(k);
         r9*cos(phi2(k))*dphi2(k)];
     
    x = Bm\Am;
    
    % save results
    dphi3(k) = x(1);
    dphi4(k) = x(2);
    dphi5(k) = x(3);
    dphi6(k) = x(4);
    dphi7(k) = x(5);
    dphi8(k) = x(6);
    
    
    
    % *** acceleration analysis ***
    
   Bm = [-r3*sin(phi3(k)),  r4*sin(phi4(k)),0,0,0,0;
         r3*cos(phi3(k)), -r4*cos(phi4(k)),0,0,0,0;
         0,0,-r5*sin(phi5(k)), -r6*sin(phi6(k)), r7*sin(phi7(k)),-r8*sin(phi8(k));
         0,0,r5*cos(phi5(k)), r6*cos(phi6(k)), -r7*cos(phi7(k)),r8*cos(phi8(k));
         -r3*sin(phi3(k)),0,0, r6*sin(phi6(k)), r10*sin(phi7(k)),0;
         r3*cos(phi3(k)),0,0, -r6*cos(phi6(k)), -r10*cos(phi7(k)),0];
    Am = [r2*cos(phi2(k))*dphi2(k)^2+r2*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r4*cos(phi4(k))*dphi4(k)^2;
         r2*sin(phi2(k))*dphi2(k)^2-r2*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r4*sin(phi4(k))*dphi4(k)^2;
         0;
         0;
         -r9*cos(phi2(k))*dphi2(k)^2-r9*sin(phi2(k))*ddphi2(k)+r3*cos(phi3(k))*dphi3(k)^2-r6*cos(phi6(k))*dphi6(k)^2-r10*cos(phi7(k))*dphi7(k)^2;
         -r9*sin(phi2(k))*dphi2(k)^2+r9*cos(phi2(k))*ddphi2(k)+r3*sin(phi3(k))*dphi3(k)^2-r6*sin(phi6(k))*dphi6(k)^2-r10*sin(phi7(k))*dphi7(k)^2];
    
    x = Bm\Am;
    % save results
    ddphi3(k) = x(1);
    ddphi4(k) = x(2);
    ddphi5(k) = x(3);
    ddphi6(k) = x(4);
    ddphi7(k) = x(5);
    ddphi8(k) = x(6);
    
    
    % *** calculate initial values for next iteration step ***
    phi3_init = phi3(k)+Ts*dphi3(k);
    phi4_init = phi4(k)+Ts*dphi4(k);
    
    
end % loop over positions



% *** create movie ***

% point P = fixed
B = 0;
% point S = fixed
A = r1*exp(j*phi1);
% define which positions we want as frames in our movie
frames = 40;    % number of frames in movie
delta = floor(t_size/frames); % time between frames
index_vec = [1:delta:t_size]';

% Create a window large enough for the whole mechanisme in all positions, to prevent scrolling.
% This is done by plotting a diagonal from (x_left, y_bottom) to (x_right, y_top), setting the
% axes equal and saving the axes into "movie_axes", so that "movie_axes" can be used for further
% plots.
x_left = -0.5*r2;
y_bottom = -0.5*max(r2,r4);
x_right = 1*r1;  %r1+1.5*r4;
y_top =  1.5*r1;  %1.5*max(r2,r4);

figure(10)
hold on
plot([x_left, x_right], [y_bottom, y_top]);
axis equal;
movie_axes = axis;   %save current axes into movie_axes

% draw and save movie frame
for m=1:length(index_vec)
    index = index_vec(m);
    
    
    C = A + r2 * exp(j*phi2(index));
    E1 = C + r3 * exp(j*phi3(index));
    E2 = B + r4 * exp(j*phi4(index));
    B2 = A - r1 * exp(j*pi/2);
    
    loop1 = [C E1 E2 B B2 A C];
    
    F1 = E1 - r10 * exp(j*phi7(index));
    D = C + r9 * exp(j*phi2(index));
    F2 = D + r6 * exp(j*phi6(index));
    
    loop3 = [C E1 F1 F2 D C];
    
    G = D - r5 * exp(j*phi5(index));
    H1 = G - r8 * exp(j*phi8(index));
    H2 = F1 - r7*exp(j*phi7(index));
    
    loop2 = [D F1 H1 H2 G D];
    
    figure(10)
    clf
    hold on
    plot(real(loop1),imag(loop1),'-o')
    plot(real(loop2),imag(loop2),'-o')
    plot(real(loop3),imag(loop3),'-o')
    
    axis(movie_axes);     % set axes as in movie_axes
    Movie(m) = getframe;  % save frame to a variable Film
end

% save movie
save fourbar_movie Movie
%close(10)


% *** plot figures ***

if fig_kin_4bar
    
    %plot assembly at a certain timestep 
    index = 1; %select 1st timestep
    B = 0;
    A = r1*exp(j*phi1);
    C = A + r2 * exp(j*phi2(index));
    E1 = C + r3 * exp(j*phi3(index));
    E2 = B + r4 * exp(j*phi4(index));
    B2 = A - r1 * exp(j*pi/2);
    F1 = E1 - r10 * exp(j*phi7(index));
    D = C + r9 * exp(j*phi2(index));
    F2 = D + r6 * exp(j*phi6(index));
    G = D - r5 * exp(j*phi5(index));
    H1 = G - r8 * exp(j*phi8(index));
    H2 = F1 - r7*exp(j*phi7(index));
    
    figure
    assembly=[A C E1 E2 F1 F2 D C E1 E2 F1 F2 D G H1 H2 E1 E2 B B2 A] ;
    plot(real(assembly),imag(assembly),'ro-')
    xlabel('[m]')
    ylabel('[m]')
    title('assembly')
    axis equal
    
    figure
    subplot(711)
    plot(t,phi2)
    ylabel('\phi_2 [rad]')
    subplot(712)
    plot(t,phi3)
    ylabel('\phi_3 [rad]')
    subplot(713)
    plot(t,phi4)
    ylabel('\phi_4 [rad]')
    subplot(714)
    plot(t,phi5)
    ylabel('\phi_5 [rad]')
    subplot(715)
    plot(t,phi6)
    ylabel('\phi_6 [rad]')
    subplot(716)
    plot(t,phi7)
    ylabel('\phi_7 [rad]')
    %subplot(717)
    figure
    plot(t,phi8)
    ylabel('\phi_8 [rad]')
    
    xlabel('t [s]')
    
    figure
    subplot(711)
    plot(t,dphi2)
    ylabel('d\phi_2 [rad/s]')
    subplot(712)
    plot(t,dphi3)
    ylabel('d\phi_3 [rad/s]')
    subplot(713)
    plot(t,dphi4)
    ylabel('d\phi_4 [rad/s]')
    subplot(714)
    plot(t,dphi5)
    ylabel('d\phi_5 [rad/s]')
    subplot(715)
    plot(t,dphi6)
    ylabel('d\phi_6 [rad/s]')
    subplot(716)
    plot(t,dphi7)
    ylabel('d\phi_7 [rad/s]')
    subplot(717)
    plot(t,dphi8)
    ylabel('d\phi_8 [rad/s]')
    xlabel('t [s]')
    
    figure
    subplot(711)
    plot(t,ddphi2)
    ylabel('dd\phi_2 [rad/s^2]')
    subplot(712)
    plot(t,ddphi3)
    ylabel('dd\phi_3 [rad/s^2]')
    subplot(713)
    plot(t,ddphi4)
    ylabel('dd\phi_4 [rad/s^2]')
    subplot(714)
    plot(t,ddphi5)
    ylabel('dd\phi_5 [rad/s^2]')
    subplot(715)
    plot(t,ddphi6)
    ylabel('dd\phi_6 [rad/s^2]')
    subplot(716)
    plot(t,ddphi7)
    ylabel('dd\phi_7 [rad/s^2]')
    subplot(717)
    plot(t,ddphi8)
    ylabel('dd\phi_8 [rad/s^2]')
    
    xlabel('t [s]')
end



