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


function [F_A_x,F_B_x,F_C_x,F_D_x_2,F_D_x_5,F_E_x_3,F_E_x_4,F_F_x,F_G_x,F_H_x,F_A_y,F_B_y,F_C_y,F_D_y_2,F_D_y_5,F_E_y_3,F_E_y_4,F_F_y,F_G_y,F_H_y,M_A] = ...
dynamics_4bar(phi2,phi3,phi4,phi5,phi6,phi7,phi8,dphi2,dphi3,dphi4,dphi5,dphi6,dphi7,dphi8,ddphi2,ddphi3,ddphi4,ddphi5,ddphi6,ddphi7,ddphi8,r2,r3,r4,r5,r6,r7,r8,r9,r10, m2,m3,m4,m5,m6,m7,m8,X2,X3,X4,X5,X6,X7,X8,Y2,Y3,Y4,Y5,Y6,Y7,Y8,J2,J3,J4,J5,J6,J7,J8,Jcog2,Jcog3,Jcog4,Jcog5,Jcog6,Jcog7,Jcog8,t,fig_dyn_4bar)


% a lot of definitions to make the matrix A and B a bit clear.
% skip the definitions for now (move down to "force analysis")
% and check them when you need them.


% cogi_P_x, cogn_P_y = vector from the centre of gravity of bar i to point P
%opm: alle Yi = 0
cog2_A_x = -X2*cos(phi2)-Y2*cos(phi2+pi/2);
cog2_A_y = -X2*sin(phi2)-Y2*sin(phi2+pi/2);
cog2_D_x = X2*cos(phi2)-Y2*cos(phi2+pi/2);
cog2_D_y = X2*sin(phi2)-Y2*sin(phi2+pi/2);
cog2_C_x = -(X2-r2)*cos(phi2)-Y2*cos(phi2+pi/2);
cog2_C_y = -(X2-r2)*sin(phi2)-Y2*sin(phi2+pi/2);
cog3_C_x = -X3*cos(phi3)-Y3*cos(phi3+pi/2);
cog3_C_y = -X3*sin(phi3)-Y3*sin(phi3+pi/2);
cog3_E_x = X3*cos(phi3)-Y3*cos(phi3+pi/2);
cog3_E_y = X3*sin(phi3)-Y3*sin(phi3+pi/2);
cog4_B_x = -X4*cos(phi4)-Y4*cos(phi4+pi/2);
cog4_B_y = -X4*sin(phi4)-Y4*sin(phi4+pi/2);
cog4_E_x = X4*cos(phi4)-Y4*cos(phi4+pi/2);
cog4_E_y = X4*sin(phi4)-Y4*sin(phi4+pi/2);
cog5_G_x = -X5*cos(phi5)-Y5*cos(phi5+pi/2);
cog5_G_y = -X5*sin(phi5)-Y5*sin(phi5+pi/2);
cog5_D_x = X5*cos(phi5)-Y5*cos(phi5+pi/2);
cog5_D_y = X5*sin(phi5)-Y5*sin(phi5+pi/2);
cog6_F_x = X6*cos(phi6)-Y6*cos(phi6+pi/2);
cog6_F_y = X6*sin(phi6)-Y6*sin(phi6+pi/2);
cog6_D_x = -X6*cos(phi6)-Y6*cos(phi6+pi/2);
cog6_D_y = -X6*sin(phi6)-Y6*sin(phi6+pi/2);
cog7_H_x = -X7*cos(phi7)-Y7*cos(phi7+pi/2);
cog7_H_y = -X7*sin(phi7)-Y7*sin(phi7+pi/2);
cog7_E_x = X7*cos(phi7)-Y7*cos(phi7+pi/2);
cog7_E_y = X7*sin(phi7)-Y7*sin(phi7+pi/2);
cog7_F_x = -(X7-r7)*cos(phi7)-Y7*cos(phi7+pi/2);
cog7_F_y = -(X7-r7)*sin(phi7)-Y7*sin(phi7+pi/2);
cog8_G_x = X8*cos(phi8)-Y8*cos(phi8+pi/2);
cog8_G_y = X8*sin(phi8)-Y8*sin(phi8+pi/2);
cog8_H_x = -X8*cos(phi8)-Y8*cos(phi8+pi/2);
cog8_H_y = -X8*sin(phi8)-Y8*sin(phi8+pi/2);

% 3D omega (dphi) and alpha (ddphi) vectors)
omega2 = [zeros(size(phi2)) zeros(size(phi2)) dphi2];
omega3 = [zeros(size(phi2)) zeros(size(phi2)) dphi3];
omega4 = [zeros(size(phi2)) zeros(size(phi2)) dphi4];
omega5 = [zeros(size(phi2)) zeros(size(phi2)) dphi5];
omega6 = [zeros(size(phi2)) zeros(size(phi2)) dphi6];
omega7 = [zeros(size(phi2)) zeros(size(phi2)) dphi7];
omega8 = [zeros(size(phi2)) zeros(size(phi2)) dphi8];
alpha2 = [zeros(size(phi2)) zeros(size(phi2)) ddphi2];
alpha3 = [zeros(size(phi2)) zeros(size(phi2)) ddphi3];
alpha4 = [zeros(size(phi2)) zeros(size(phi2)) ddphi4];
alpha5 = [zeros(size(phi2)) zeros(size(phi2)) ddphi5];
alpha6 = [zeros(size(phi2)) zeros(size(phi2)) ddphi6];
alpha7 = [zeros(size(phi2)) zeros(size(phi2)) ddphi7];
alpha8 = [zeros(size(phi2)) zeros(size(phi2)) ddphi8];

% 3D model vectors

A_cog2_vec = [-cog2_A_x    -cog2_A_y    zeros(size(phi2))];
% D_cog2_vec = [-cog2_D_x    -cog2_D_y    zeros(size(phi2))];
% C_cog2_vec = [-cog2_C_x    -cog2_C_y    zeros(size(phi2))];
C_cog3_vec = [-cog3_C_x    -cog3_C_y    zeros(size(phi2))];
B_cog4_vec = [-cog4_B_x    -cog4_B_y    zeros(size(phi2))];
D_cog6_vec = [-cog6_D_x    -cog6_D_y    zeros(size(phi2))];
H_cog8_vec = [-cog8_H_x    -cog8_H_y    zeros(size(phi2))];
G_cog5_vec = [-cog5_G_x    -cog5_G_y    zeros(size(phi2))];
H_cog7_vec = [-cog7_H_x    -cog7_H_y    zeros(size(phi2))];
AD_vec = [(r2+r9)*cos(phi2) (r2+r9)*sin(phi2) zeros(size(phi2))];
AC_vec = [r2*cos(phi2) r2*sin(phi2) zeros(size(phi2))];
CE_vec = [r3*cos(phi3) r3*sin(phi3) zeros(size(phi2))];
DF_vec = [r6*cos(phi6) r6*sin(phi6) zeros(size(phi2))];
FH_vec = [-r7*cos(phi7) -r7*sin(phi7) zeros(size(phi2))];
HG_vec = [r8*cos(phi8) r8*sin(phi8) zeros(size(phi2))];


% acceleration vectors
acc_2 =       cross(omega2,cross(omega2,A_cog2_vec))+cross(alpha2,A_cog2_vec);
acc_D =       cross(omega2,cross(omega2,AD_vec    ))+cross(alpha2,AD_vec    );
acc_C =       cross(omega2,cross(omega2,AC_vec    ))+cross(alpha2,AC_vec    );
acc_3 = acc_C+cross(omega3,cross(omega3,C_cog3_vec))+cross(alpha3,C_cog3_vec);
acc_E = acc_C+cross(omega3,cross(omega3,CE_vec    ))+cross(alpha3,CE_vec);
acc_4 =       cross(omega4,cross(omega4,B_cog4_vec))+cross(alpha4,B_cog4_vec);


acc_6 = acc_D+cross(omega6,cross(omega6,D_cog6_vec))+cross(alpha6,D_cog6_vec);
acc_F = acc_D+cross(omega6,cross(omega6,DF_vec))+cross(alpha6,DF_vec);
acc_H = acc_F+cross(omega7,cross(omega7,FH_vec))+cross(alpha7,FH_vec);
acc_8 = acc_H+cross(omega8,cross(omega8,H_cog8_vec))+cross(alpha8,H_cog8_vec);
acc_G = acc_H+cross(omega8,cross(omega8,HG_vec))+cross(alpha8,HG_vec);
acc_5 = acc_G+cross(omega5,cross(omega5,G_cog5_vec))+cross(alpha5,G_cog5_vec);
acc_7 = acc_H+cross(omega7,cross(omega7,H_cog7_vec))+cross(alpha7,H_cog7_vec);

acc_2x = acc_2(:,1);
acc_2y = acc_2(:,2);
acc_3x = acc_3(:,1);
acc_3y = acc_3(:,2);
acc_4x = acc_4(:,1);
acc_4y = acc_4(:,2);
acc_5x = acc_5(:,1);
acc_5y = acc_5(:,2);
acc_6x = acc_6(:,1);
acc_6y = acc_6(:,2);
acc_7x = acc_7(:,1);
acc_7y = acc_7(:,2);
acc_8x = acc_8(:,1);
acc_8y = acc_8(:,2);

%Velocity vectors
vel_2 = cross(omega2, A_cog2_vec);
vel_C = cross(omega2,AC_vec);
vel_D = cross(omega2, AD_vec);
vel_4 = cross(omega4, B_cog4_vec);
vel_3 = vel_C + cross(omega3, C_cog3_vec);
vel_6 = vel_D + cross(omega6, D_cog6_vec);
vel_F = vel_D + cross(omega6, DF_vec);
vel_H = vel_F + cross(omega7, FH_vec);
vel_7 = vel_H + cross(omega7, H_cog7_vec);
vel_8 = vel_H + cross(omega8, H_cog8_vec);
vel_G = vel_H + cross(omega8, HG_vec);
vel_5 = vel_G + cross(omega5, G_cog5_vec);

vel_2x = vel_2(:,1);
vel_2y = vel_2(:,2);
vel_3x = vel_3(:,1);
vel_3y = vel_3(:,2);
vel_4x = vel_4(:,1);
vel_4y = vel_4(:,2);
vel_5x = vel_5(:,1);
vel_5y = vel_5(:,2);
vel_6x = vel_6(:,1);
vel_6y = vel_6(:,2);
vel_7x = vel_7(:,1);
vel_7y = vel_7(:,2);
vel_8x = vel_8(:,1);
vel_8y = vel_8(:,2);





% **********************
% *** force analysis ***
% **********************

% allocate matrices for force (F) and moment (M)
    F_A_x = zeros(size(phi2));
    F_A_y = zeros(size(phi2));
    F_B_x = zeros(size(phi2));
    F_B_y = zeros(size(phi2));
    F_C_x =zeros(size(phi2));
    F_C_y = zeros(size(phi2));
    F_D_x_2 =zeros(size(phi2));
    F_D_y_2 = zeros(size(phi2));
    F_D_x_5 = zeros(size(phi2));
    F_D_y_5 = zeros(size(phi2));
    F_E_x_3 = zeros(size(phi2));
    F_E_y_3 = zeros(size(phi2));
    F_E_x_4 = zeros(size(phi2));
    F_E_y_4 = zeros(size(phi2));
    F_F_x = zeros(size(phi2));
    F_F_y = zeros(size(phi2));
    F_G_x = zeros(size(phi2));
    F_G_y = zeros(size(phi2));
    F_H_x =zeros(size(phi2));
    F_H_y = zeros(size(phi2));
    M_A  = zeros(size(phi2));


% calculate dynamics for each time step
t_size = size(t,1);    % number of simulation steps
for k=1:t_size

  
    
    
Am = [1, 0,	0,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
    0,	1,	0,	0,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0;
    -cog2_A_y(k),	cog2_A_x(k),	0,	0,	-cog2_C_y(k),	cog2_C_x(k),	-cog2_D_y(k),	cog2_D_x(k),	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1;
    0	0	0	0	-1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	cog3_C_y(k)	-cog3_C_x(k)	0	0	0	0	-cog3_E_y(k)	cog3_E_x(k)	0	0	0	0	0	0	0	0	0;
    0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0;
    0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0;
    0	0	-cog4_B_y(k)	cog4_B_x(k)	0	0	0	0	0	0	0 0 -cog4_E_y(k)	cog4_E_x(k)		0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0;
    0	0	0	0	0	0	0	0	-cog5_D_y(k)	cog5_D_x(k)	0	0	0	0	0	0	-cog5_G_y(k)	cog5_G_x(k)	0	0	0;
    0	0	0	0	0	0	-1	0	-1	0	0	0	0	0	1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	-1	0	-1	0	0	0	0	0	1	0	0	0	0	0;
    0	0	0	0	0	0	cog6_D_y(k)	-cog6_D_x(k)	cog6_D_y(k)	-cog6_D_x(k)	0	0	0	0	-cog6_F_y(k)	cog6_F_x(k)	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	0	-1	0	-1	0	0	0	1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	-1	0	-1	0	-1	0	0	0	1	0;
    0	0	0	0	0	0	0	0	0	0	cog7_E_y(k)	-cog7_E_x(k)	cog7_E_y(k)	-cog7_E_x(k)	cog7_F_y(k)	-cog7_F_x(k)	0	0	-cog7_H_y(k)	cog7_H_x(k)	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	-1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	-1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	cog8_G_y(k)	-cog8_G_x(k)	cog8_H_y(k)	-cog8_H_x(k)	0];

    
    
    
    
    
    
    
    
    
    
    
  Bm = [ m2*acc_2x(k);
        m2*acc_2y(k);
        Jcog2*ddphi2(k);
        m3*acc_3x(k);
        m3*acc_3y(k);
        Jcog3*ddphi3(k);
        m4*acc_4x(k);
        m4*acc_4y(k);
        Jcog4*ddphi4(k);
        m5*acc_5x(k);
        m5*acc_5y(k);
        Jcog5*ddphi5(k);
        m6*acc_6x(k);
        m6*acc_6y(k);
        Jcog6*ddphi6(k);
        m7*acc_7x(k);
        m7*acc_7y(k);
        Jcog7*ddphi7(k);
        m8*acc_8x(k);
        m8*acc_8y(k);
        Jcog8*ddphi8(k)];
    Am;
   [R, bj] = rref(Am);
    x = Am\Bm;
    
    % save results
    F_A_x(k) = x(1);
    F_A_y(k) = x(2);
    F_B_x(k) = x(3);
    F_B_y(k) = x(4);
    F_C_x(k) = x(5);
    F_C_y(k) = x(6);
    F_D_x_2(k) = x(7);
    F_D_y_2(k) = x(8);
    F_D_x_5(k) = x(9);
    F_D_y_5(k) = x(10);
    F_E_x_3(k) = x(11);
    F_E_y_3(k) = x(12);
    F_E_x_4(k) = x(13);
    F_E_y_4(k) = x(14);
    F_F_x(k) = x(15);
    F_F_y(k) = x(16);
    F_G_x(k) = x(17);
    F_G_y(k) = x(18);
    F_H_x(k) = x(19);
    F_H_y(k) = x(20);
    M_A(k)   = x(21);
end



% **********************
% *** plot figures ***
% **********************

if fig_dyn_4bar
    
    figure
    subplot(221)
    plot(F_A_x,F_A_y),grid
    xlabel('F_A_x [N]')
    ylabel('F_A_y [N]')
    axis tight
    subplot(222)
    plot(F_B_x,F_B_y),grid
    xlabel('F_Q_x [N]')
    ylabel('F_Q_y [N]')
    axis tight
    subplot(223)
    plot(F_C_x,F_C_y),grid
    xlabel('F_R_x [N]')
    ylabel('F_R_y [N]')
    axis tight
    
    
    figure
    subplot(221)
    plot(F_D_x_2,F_D_y_2),grid
    xlabel('F_D_x_2 [N]')
    ylabel('F_D_y_2 [N]')
    axis tight
    subplot(222)
    plot(F_D_x_5,F_D_y_5),grid
    xlabel('F_D_x_5 [N]')
    ylabel('F_D_y_5 [N]')
    axis tight
    subplot(223)
    plot((-F_D_x_2 -F_D_x_5) ,(-F_D_y_2-F_D_y_5)),grid
    xlabel('F_D_x_6 [N]')
    ylabel('F_D_y_6 [N]')
    axis tight
    
    
    
    figure
    subplot(221)
    plot(F_E_x_3,F_E_y_3),grid
    xlabel('F_E_x_3 [N]')
    ylabel('F_E_y_3 [N]')
    axis tight
    subplot(222)
    plot(F_E_x_4,F_E_y_4),grid
    xlabel('F_E_x_4 [N]')
    ylabel('F_E_y_4 [N]')
    axis tight
    subplot(223)
    plot(-F_E_x_3 -F_E_x_4 ,-F_E_y_3-F_E_y_4),grid
    xlabel('F_E_x_7 [N]')
    ylabel('F_E_y_7 [N]')
    axis tight
    
    
    
    figure
    subplot(221)
    plot(F_F_x,F_F_y),grid
    xlabel('F_F_x [N]')
    ylabel('F_F_y [N]')
    axis tight
    subplot(222)
    plot(F_G_x,F_G_y),grid
    xlabel('F_G_x [N]')
    ylabel('F_G_y [N]')
    axis tight
    subplot(223)
    plot(F_H_x,F_H_y),grid
    xlabel('F_H_x [N]')
    ylabel('F_H_y [N]')
    axis tight

    
    figure
    plot(t,M_A)
    ylabel('M_A [N-m]')
    xlabel('t [s]')
    
end

% *************************
% *** controle dynamica ***
% *************************

M_controle = zeros(t_size);

dEdt = zeros(t_size);
for k=1:t_size

    dEdt(k) = Jcog2*dphi2(k)*ddphi2(k) + Jcog3*dphi3(k)*ddphi3(k) + Jcog4*dphi4(k)*ddphi4(k) + ...
        Jcog5*dphi5(k)*ddphi5(k) + Jcog6*dphi6(k)*ddphi6(k) + Jcog7*dphi7(k)*ddphi7(k) + Jcog8*dphi8(k)*ddphi8(k)+ ...
        m2*(vel_2x(k)*acc_2x(k)+vel_2y(k)*acc_2y(k))+m3*(vel_3x(k)*acc_3x(k)+vel_3y(k)*acc_3y(k))+m4*(vel_4x(k)*acc_4x(k)+vel_4y(k)*acc_4y(k))+m5*(vel_5x(k)*acc_5x(k)+vel_5y(k)*acc_5y(k))+...
        m6*(vel_6x(k)*acc_6x(k)+vel_6y(k)*acc_6y(k))+m7*(vel_7x(k)*acc_7x(k)+vel_7y(k)*acc_7y(k))+m8*(vel_8x(k)*acc_8x(k)+vel_8y(k)*acc_8y(k));
    M_controle(k) = dEdt(k)/dphi2(k);
end

figure
subplot(211)
plot(t,M_controle(:,1))
ylabel('M_{controle} [N-m]')
xlabel('t [s]')
subplot(212)
plot(t,M_A)
ylabel('M_A [N-m]')
xlabel('t [s]')

figure
plot(t, M_A - M_controle(:,1))
ylabel('M_A - M_{controle} [N-m]')
xlabel('t [s]')






