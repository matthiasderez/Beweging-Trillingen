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

clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data initialization (all data is converted to SI units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program data
fig_kin_4bar = 1;        % draw figures of kinematic analysis if 1
fig_dyn_4bar = 1;        % draw figures of dynamic analysis if 1

% kinematic parameters (link lengths)
r1 = 1; %BA
r2 = 0.4; %AC
r3 = 0.4; %CE
r4 = 1; %BE
r5 = 1; %GD
r6 = 0.4; %DF
r7 = 0.4; %HF
r8 = 1; %HG
r9 = 0.6; %DC
r10 = 0.6; %FE
phi1 = pi/2; %hoek horizontale met AB

% dynamic parameters, defined in a local frame on each of the bars.
X2 = (r2+r9)/2;               % X coordinates of cog (centre of gravity)
X3 = r3/2;
X4 = r4/2;
X5 = r5/2;
X6 = r6/2;
X7 = (r7+r10)/2;
X8 = r8/2;

Y2 = 0;                  % Y coordinates of cog
Y3 = 0;
Y4 = 0;
Y5 = 0;
Y6 = 0;
Y7 = 0;
Y8 = 0;





m2 = (r2+r9)*3.12; %staal met doorsnede 4cm² --> 3.12kg/m
m3 = r3*3.12;
m4 = r4*3.12;
m5 = r5*3.12;
m6 = r6*3.12;
m7 = (r7+r10)*3.12;
m8 = r8*3.12;

J2 = m2*(r2+r9)^2/12 + m2*( X2^2 + Y2^2);
J3 = m3*r3^2/12 + m3*( X3^2 + Y3^2);
J4 = m4*r4^2/12+ m4*( X4^2 + Y4^2);
J5 = m5*r5^2/12+ m5*( X5^2 + Y5^2);
J6 = m6*r6^2/12+ m6*( X6^2 + Y6^2);
J7 = m7*(r7+r10)^2/12+ m7*( X7^2 + Y7^2);
J8 = m8*r8^2/12+ m8*( X8^2 + Y8^2);


Jcog2 = m2*(r2+r9)^2/12;
Jcog3 = m3*r3^2/12;
Jcog4 = m4*r4^2/12;
Jcog5 = m5*r5^2/12;
Jcog6 = m6*r6^2/12;
Jcog7 = m7*(r7+r10)^2/12;
Jcog8 = m8*r8^2/12;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1. Determination of Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t_begin = 0;                   % start time of simulation
t_end = 10;                    % end time of simulation
Ts = 0.05;                     % time step of simulation
t = [t_begin:Ts:t_end]';       % time vector

% initialization of driver
omega = 0.5;
A = pi/9;
phi2=3*pi/2+A+A/2*sin(omega*t);
dphi2=omega*A/2*sin(omega*t);
ddphi2=omega^2*A/2*cos(omega*t);

% position analysis
epsilon = 0.05;
phi3_init = A;    % initial condition for first step of position analysis with fsolve (phi3 and phi4)
phi4_init = pi/2 - A;  % VERY IMPORTANT because it determines which branch of the mechanism you're in
phi5_init = 3*pi/2 - A;
phi6_init = A;
phi7_init = pi/2 + A;
phi8_init = pi/2-epsilon;


% calculation of the kinematics (see kin_4bar.m)
[phi3,phi4,phi5,phi6,phi7,phi8,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8] = kinematics_4bar(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,phi1,phi2,dphi2,ddphi2,phi3_init,phi4_init,phi5_init,phi6_init,phi7_init,phi8_init,t,fig_kin_4bar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2. Dynamics Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculation of the dynamics (see dyn_4bar.m)
[F_A_x,F_B_x,F_C_x,F_D_x_2,F_D_x_5,F_E_x_3,F_E_x_4,F_F_x,F_G_x,F_H_x,F_A_y,F_B_y,F_C_y,F_D_y_2,F_D_y_5,F_E_y_3,F_E_y_4,F_F_y,F_G_y,F_H_y,M_A] = dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,r2,r3,r4,r5,r6,r7,r8,r9,r10,m2,m3,m4,m5,m6,m7,m8,X2,X3,X4,X5,X6,X7,X8,Y2,Y3,Y4,Y5,Y6,Y7,Y8,J2,J3,J4,J5,J6,J7,J8,Jcog2,Jcog3,Jcog4,Jcog5,Jcog6,Jcog7,Jcog8,t,fig_dyn_4bar);
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3. Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
load fourbar_movie Movie
movie(Movie)

