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


function F=loop_closure_eqs(phi_init,phi2,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,phi1)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi3=phi_init(1);
phi4=phi_init(2);
phi5=phi_init(3);
phi6=phi_init(4);
phi7=phi_init(5);
phi8=phi_init(6);


% loop closure equations:
F(1)=r1*cos(phi1)+r2*cos(phi2)+r3*cos(phi3)-r4*cos(phi4);
F(2)=r1*sin(phi1)+r2*sin(phi2)+r3*sin(phi3)-r4*sin(phi4);
F(3)=r8*cos(phi8)+r5*cos(phi5)+r6*cos(phi6)-r7*cos(phi7);
F(4)=r8*sin(phi8)+r5*sin(phi5)+r6*sin(phi6)-r7*sin(phi7);
F(5)=r3*cos(phi3)-r10*cos(phi7)-r6*cos(phi6)-r9*cos(phi2);
F(6)=r3*sin(phi3)-r10*sin(phi7)-r6*sin(phi6)-r9*sin(phi2);

