clear;close all;clc
%Denavit Hartenberg for Guitar Stand

syms a_0 alpha_0
syms d_1 theta_1 a_1 alpha_1 
syms d_2 theta_2 a_2 alpha_2
syms d_3 theta_3 a_3 alpha_3
syms d_4 theta_4 a_4 alpha_4
syms d_5 theta_5

%/////////////////////////////////////////////////////////////////////////////////////////////////////
%Establishing D-H Matrix  

T_1_0_row_1=[cos(theta_1) -sin(theta_1) 0 a_0];
T_1_0_row_2=[sin(theta_1)*cos(alpha_0) cos(theta_1)*cos(alpha_0) -sin(alpha_0) -d_1*sin(alpha_0)];
T_1_0_row_3=[sin(theta_1)*sin(alpha_0) cos(theta_1)*sin(alpha_0) cos(alpha_0) d_1*cos(alpha_0)];
T_1_0_row_4=[0 0 0 1];
T_1_0=[T_1_0_row_1;T_1_0_row_2;T_1_0_row_3;T_1_0_row_4];

%Establishing D-H Matrix 2
T_2_1_row_1=[cos(theta_2) -sin(theta_2) 0 a_1];
T_2_1_row_2=[sin(theta_2)*cos(alpha_1) cos(theta_2)*cos(alpha_1) -sin(alpha_1) -d_2*sin(alpha_1)];
T_2_1_row_3=[sin(theta_2)*sin(alpha_1) cos(theta_2)*sin(alpha_1) cos(alpha_1) d_2*cos(alpha_1)];
T_2_1_row_4=[0 0 0 1];
T_2_1=[T_2_1_row_1;T_2_1_row_2;T_2_1_row_3;T_2_1_row_4];


%Establishing D-H Matrix 3

T_3_2_row_1=[cos(theta_3) -sin(theta_3) 0 a_2];
T_3_2_row_2=[sin(theta_3)*cos(alpha_2) cos(theta_3)*cos(alpha_2) -sin(alpha_2) -d_3*sin(alpha_2)];
T_3_2_row_3=[sin(theta_3)*sin(alpha_2) cos(theta_3)*sin(alpha_2) cos(alpha_2) d_3*cos(alpha_2)];
T_3_2_row_4=[0 0 0 1];
T_3_2=[T_3_2_row_1;T_3_2_row_2;T_3_2_row_3;T_3_2_row_4];

%Establishing D-H Matrix 4

T_4_3_row_1=[cos(theta_4) -sin(theta_4) 0 a_3];
T_4_3_row_2=[sin(theta_4)*cos(alpha_3) cos(theta_4)*cos(alpha_3) -sin(alpha_3) -d_4*sin(alpha_3)];
T_4_3_row_3=[sin(theta_4)*sin(alpha_3) cos(theta_4)*sin(alpha_3) cos(alpha_3) d_4*cos(alpha_3)];
T_4_3_row_4=[0 0 0 1];
T_4_3=[T_4_3_row_1;T_4_3_row_2;T_4_3_row_3;T_4_3_row_4];

%Establishing D-H Matrix 5

T_5_4_row_1=[cos(theta_5) -sin(theta_5) 0 a_4];
T_5_4_row_2=[sin(theta_5)*cos(alpha_4) cos(theta_5)*cos(alpha_4) -sin(alpha_4) -d_5*sin(alpha_4)];
T_5_4_row_3=[sin(theta_5)*sin(alpha_4) cos(theta_5)*sin(alpha_4) cos(alpha_4) d_5*cos(alpha_4)];
T_5_4_row_4=[0 0 0 1];
T_5_4=[T_5_4_row_1;T_5_4_row_2;T_5_4_row_3;T_5_4_row_4];

%Establishing T Matrices

T_2_0=T_1_0*T_2_1;
T_3_0=T_1_0*T_2_1*T_3_2;
T_4_0=T_1_0*T_2_1*T_3_2*T_4_3;
T_5_0=T_1_0*T_2_1*T_3_2*T_4_3*T_5_4;


%/////////////////////////////////////////////////////////////////////////////////////////////////////////////
%Creating Jacobian Matrices

q1=d_1;
q2=d_2;
q3=theta_3;
q4=theta_4;
q=[q1 q2 q3 q4];

r_x=T_5_0(1,4);
r_y=T_5_0(2,4);
r_z=T_5_0(3,4);
T_homogeneous=T_5_0(4,:);

J_drx_dq1=diff(r_x,q1);
J_drx_dq2=diff(r_x,q2);
J_drx_dq3=diff(r_x,q3);
J_drx_dq4=diff(r_x,q4);

J_dry_dq1=diff(r_y,q1);
J_dry_dq2=diff(r_y,q2);
J_dry_dq3=diff(r_y,q3);
J_dry_dq4=diff(r_y,q4);

J_drz_dq1=diff(r_z,q1);
J_drz_dq2=diff(r_z,q2);
J_drz_dq3=diff(r_z,q3);
J_drz_dq4=diff(r_z,q4);

J_velocity_vector=[J_drx_dq1 J_drx_dq2 J_drx_dq3 J_drx_dq4; J_dry_dq1 J_dry_dq2 J_dry_dq3 J_dry_dq4 ; J_drz_dq1 J_drz_dq2 J_drz_dq3 J_drz_dq4];

epsi_1=1;
epsi_1_bar=~epsi_1;
epsi_2=1;
epsi_2_bar=~epsi_2;
epsi_3=0;
epsi_3_bar=~epsi_3;
epsi_4=0;
epsi_4_bar=~epsi_4;

z_1_1_hat=[0;0;1;1];
z_2_2_hat=[0;0;1;1];
z_3_3_hat=[0;0;1;1];
z_4_4_hat=[0;0;1;1];

z_1_0_hat=T_1_0*z_1_1_hat;
z_2_0_hat=T_2_0*z_1_1_hat;
z_3_0_hat=T_3_0*z_3_3_hat;
z_4_0_hat=T_4_0*z_4_4_hat;
omega_1=epsi_1_bar*z_1_0_hat;
omega_2=epsi_2_bar*z_2_0_hat;
omega_3=epsi_3_bar*z_3_0_hat;
omega_4=epsi_4_bar*z_4_0_hat;

J_omega_vector=[omega_1 omega_2 omega_3 omega_4];
J_omega_vector=J_omega_vector((1:3),:);

%////////////////////////////////////////////////////////////////////////////////////////////////////////

%Inputting initial values

a_0_input=0;
alpha_0_input=0;

d_1_input=2;
theta_1_input=0;
a_1_input=0;
alpha_1_input=0;

d_2_input=10;
theta_2_input=0;
a_2_input=0;
alpha_2_input=0;

d_3_input=0;
theta_3_input=0;
a_3_input=0;
alpha_3_input=90;

d_4_input=10.5;
theta_4_input=90;
a_4_input=2;
alpha_4_input=0;

d_5_input=3.5;
theta_5_input=0;

%////////////////////////////////////////////////////////////////////////////////////////////////////////

%Evaluating transformation and D-H matrices

T_1_0_eval=subs(T_1_0,[a_0,alpha_0,d_1,theta_1],[a_0_input,alpha_0_input,d_1_input,theta_1_input]);
T_2_1_eval=subs(T_2_1,[a_1,alpha_1,d_2,theta_2],[a_1_input,alpha_1_input,d_2_input,theta_2_input]);
T_3_2_eval=subs(T_3_2,[a_2,alpha_2,d_3,theta_3],[a_2_input,alpha_2_input,d_3_input,theta_3_input]);
T_4_3_eval=subs(T_4_3,[a_3,alpha_3,d_4,theta_4],[a_3_input,alpha_3_input,d_4_input,theta_4_input]);
T_5_4_eval=subs(T_5_4,[a_4,alpha_4,d_5,theta_5],[a_4_input,alpha_4_input,d_5_input,theta_5_input]);

T_2_0_eval=T_1_0_eval*T_2_1_eval;
T_3_0_eval=T_1_0_eval*T_2_1_eval*T_3_2_eval;
T_4_0_eval=T_1_0_eval*T_2_1_eval*T_3_2_eval*T_4_3_eval;
T_5_0_eval=T_1_0_eval*T_2_1_eval*T_3_2_eval*T_4_3_eval*T_5_4_eval;


%//////////////////////////////////////////////////////////////////////////////////////////////////////////

%Inverse Kinematics Iteration 1

q_initial=[2.2;5;45;45];
P_1=[-7; -3.5; 22.5];
P_1_norm=norm(P_1);
P_1_angle_rot_x=acosd(P_1(1)/P_1_norm);
P_1_angle_rot_y=acosd(P_1(2)/P_1_norm);
P_1_angle_rot_z=acosd(P_1(3)/P_1_norm);
P_1_total=[P_1;P_1_angle_rot_x;P_1_angle_rot_y;P_1_angle_rot_z]
P_1_total_norm=norm(P_1_total);
count=0;

P_interp_total=[0;0;0;0;0;0];
P_interp_total_norm=norm(P_interp_total);
error=abs((P_1_total_norm-P_interp_total_norm)/P_1_total_norm);

%Uncomment to Plot Position
x1=[P_1_total(1,1);P_interp_total(1,1)];
y1=[P_1_total(2,1);P_interp_total(2,1)];
z1=[P_1_total(3,1);P_interp_total(3,1)];


plot3(x1,y1,z1,'r-o')
xlabel('x');ylabel('y');zlabel('z');
grid on
hold on

%Uncomment to Plot Angular Velocity
% theta_x_1=[P_1_total(4);P_interp_total(4)];
% theta_y_1=[P_1_total(5);P_interp_total(5)];
% theta_z_1=[P_1_total(6);P_interp_total(6)];

% plot3(theta_x_1,theta_y_1,theta_z_1)
% xlabel('theta_x');ylabel('theta_y');zlabel('theta_z');
% grid on
% hold on

    
 while error>0.01
     
    T_5_0_interp=subs(T_5_0,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,d_5,theta_5],...
    [a_0_input,alpha_0_input,q_initial(1),theta_1_input,a_1_input,alpha_1_input,q_initial(2),theta_2_input,a_2_input,alpha_2_input,d_3_input,q_initial(3),...
    a_3_input,alpha_3_input,d_4_input,q_initial(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);
    T_5_0_interp=T_5_0_interp((1:3),:);

    P_interp=T_5_0_interp(:,4);
    P_interp_norm=norm(P_interp);
    P_interp_angle_rot_x=acosd(P_interp(1)/P_interp_norm);
    P_interp_angle_rot_y=acosd(P_interp(2)/P_interp_norm);
    P_interp_angle_rot_z=acosd(P_interp(3)/P_interp_norm);
    P_interp_total=[P_interp;P_interp_angle_rot_x;P_interp_angle_rot_y;P_interp_angle_rot_z];
    P_interp_total_norm=norm(P_interp_total);
    
    
    J_velocity_vector_eval=subs(J_velocity_vector,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,...
        d_5,theta_5],[a_0_input,alpha_0_input,q_initial(1),theta_1_input,a_1_input,alpha_1_input,q_initial(2),theta_2_input,a_2_input,alpha_2_input,d_3_input...
        q_initial(3),a_3_input,alpha_3_input,d_4_input,q_initial(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);

    J_omega_vector_eval=subs(J_omega_vector,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,...
        d_5,theta_5],[a_0_input,alpha_0_input,q_initial(1),theta_1_input,a_1_input,alpha_1_input,q_initial(2),theta_2_input,a_2_input,alpha_2_input,d_3_input...
        q_initial(3),a_3_input,alpha_3_input,d_4_input,q_initial(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);

    J_eval=[J_velocity_vector_eval;J_omega_vector_eval];
    m=pinv(J_eval);
    
    P_delta=(P_1_total-P_interp_total);
    
    
    
    q_initial=q_initial+(m*(P_delta/20));
    
    if q_initial(1)>5
        J_drx_dq1=0;
        J_drx_dq2=diff(r_x,q2);
        J_drx_dq3=diff(r_x,q3);
        J_drx_dq4=diff(r_x,q4);

        J_dry_dq1=0;
        J_dry_dq2=diff(r_y,q2);
        J_dry_dq3=diff(r_y,q3);
        J_dry_dq4=diff(r_y,q4);

        J_drz_dq1=0;
        J_drz_dq2=diff(r_z,q2);
        J_drz_dq3=diff(r_z,q3);
        J_drz_dq4=diff(r_z,q4);

        J_velocity_vector=[J_drx_dq1 J_drx_dq2 J_drx_dq3 J_drx_dq4; J_dry_dq1 J_dry_dq2 J_dry_dq3 J_dry_dq4 ; J_drz_dq1 J_drz_dq2 J_drz_dq3 J_drz_dq4];
    else
    end
    
    if q_initial(2)>15
        J_drx_dq1=0;
        J_drx_dq2=0;
        J_drx_dq3=diff(r_x,q3);
        J_drx_dq4=diff(r_x,q4);

        J_dry_dq1=0;
        J_dry_dq2=0;
        J_dry_dq3=diff(r_y,q3);
        J_dry_dq4=diff(r_y,q4);

        J_drz_dq1=0;
        J_drz_dq2=0;
        J_drz_dq3=diff(r_z,q3);
        J_drz_dq4=diff(r_z,q4);

        J_velocity_vector=[J_drx_dq1 J_drx_dq2 J_drx_dq3 J_drx_dq4; J_dry_dq1 J_dry_dq2 J_dry_dq3 J_dry_dq4 ; J_drz_dq1 J_drz_dq2 J_drz_dq3 J_drz_dq4];

    else
    end
    error=abs((P_1_total_norm-P_interp_total_norm)/P_1_total_norm)
    count=count+1;
%Uncomment to Plot position 
    x1(count,1)=P_interp_total(1,1);
    y1(count,1)=P_interp_total(2,1);
    z1(count,1)=P_interp_total(3,1);

%Uncomment to Plot Angular Velocity
%       theta_x_1(count,1)=P_interp_total(4,1);
%       theta_y_1(count,1)=P_interp_total(5,1);
%       theta_z_1(count,1)=P_interp_total(6,1);
 end
 error
 
%Uncomment to Plot Position
 plot3(x1,y1,z1,'g-o')
 title('Plotted trajectory of Position for Iteration 1');

%Uncomment to Plot Angular Velocity
% plot3(theta_x_1,theta_y_1,theta_z_1,'g-o')
% title('Plotted trajectory of Angle for Iteration 1');
% 


%//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

% Inverse Kinematics Iteration 2
q_initial_2=[3;6;53;53];
P_2=[-7.3; -1.8; 24];
P_2_norm=norm(P_2);
P_2_angle_rot_x=acosd(P_2(1)/P_2_norm);
P_2_angle_rot_y=acosd(P_2(2)/P_2_norm);
P_2_angle_rot_z=acosd(P_2(3)/P_2_norm);
P_2_total=[P_2;P_2_angle_rot_x;P_2_angle_rot_y;P_2_angle_rot_z]
P_2_total_norm=norm(P_2_total);
count=0;



P_interp_2_total=[1;1;1;0;0;0];
P_interp_2_total_norm=norm(P_interp_total);
error_2=abs((P_2_total_norm-P_interp_2_total_norm)/P_2_total_norm);

% Uncomment to plot Position
x2=[P_2_total(1,1);P_interp_2_total(1,1)];
y2=[P_2_total(2,1);P_interp_2_total(2,1)];
z2=[P_2_total(3,1);P_interp_2_total(3,1)];
plot3(x2,y2,z2,'r-o')
xlabel('x');ylabel('y');zlabel('z');
grid on
hold on



% Unncomment to Plot Angular Velocity
% theta_x_2=[P_2_total(4);P_interp_2_total(4)];
% theta_y_2=[P_2_total(5);P_interp_2_total(5)];
% theta_z_2=[P_2_total(6);P_interp_2_total(6)];
% 
% plot3(theta_x_2,theta_y_2,theta_z_2)
% xlabel('theta_x');ylabel('theta_y');zlabel('theta_z');
% grid on
% hold on


    
 while error_2>0.01
     
    T_5_0_interp_2=subs(T_5_0,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,d_5,theta_5],...
    [a_0_input,alpha_0_input,q_initial_2(1),theta_1_input,a_1_input,alpha_1_input,q_initial_2(2),theta_2_input,a_2_input,alpha_2_input,d_3_input,q_initial_2(3),...
    a_3_input,alpha_3_input,d_4_input,q_initial_2(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);
    T_5_0_interp_2=T_5_0_interp_2((1:3),:);

    P_interp_2=T_5_0_interp_2(:,4);
    P_interp_2_norm=norm(P_interp_2);
    P_interp_2_angle_rot_x=acosd(P_interp_2(1)/P_interp_2_norm);
    P_interp_2_angle_rot_y=acosd(P_interp_2(2)/P_interp_2_norm);
    P_interp_2_angle_rot_z=acosd(P_interp_2(3)/P_interp_2_norm);
    P_interp_2_total=[P_interp_2;P_interp_2_angle_rot_x;P_interp_2_angle_rot_y;P_interp_2_angle_rot_z];
    P_interp_2_total_norm=norm(P_interp_2_total);
    
    
    J_velocity_vector_eval=subs(J_velocity_vector,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,...
        d_5,theta_5],[a_0_input,alpha_0_input,q_initial_2(1),theta_1_input,a_1_input,alpha_1_input,q_initial_2(2),theta_2_input,a_2_input,alpha_2_input,d_3_input...
        q_initial_2(3),a_3_input,alpha_3_input,d_4_input,q_initial_2(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);

    J_omega_vector_eval=subs(J_omega_vector,[a_0,alpha_0,d_1,theta_1,a_1,alpha_1,d_2,theta_2,a_2,alpha_2,d_3,theta_3,a_3,alpha_3,d_4,theta_4,a_4,alpha_4,...
        d_5,theta_5],[a_0_input,alpha_0_input,q_initial_2(1),theta_1_input,a_1_input,alpha_1_input,q_initial_2(2),theta_2_input,a_2_input,alpha_2_input,d_3_input...
        q_initial_2(3),a_3_input,alpha_3_input,d_4_input,q_initial_2(4),a_4_input,alpha_4_input,d_5_input,theta_5_input]);

    J_eval=[J_velocity_vector_eval;J_omega_vector_eval];
    m=pinv(J_eval);
    
    P_delta_2=(P_2_total-P_interp_2_total);
    
   
    
    q_initial_2=q_initial_2+(m*(P_delta_2/20));
    
    if q_initial_2(1)>5
        J_drx_dq1=0;
        J_drx_dq2=diff(r_x,q2);
        J_drx_dq3=diff(r_x,q3);
        J_drx_dq4=diff(r_x,q4);

        J_dry_dq1=0;
        J_dry_dq2=diff(r_y,q2);
        J_dry_dq3=diff(r_y,q3);
        J_dry_dq4=diff(r_y,q4);

        J_drz_dq1=0;
        J_drz_dq2=diff(r_z,q2);
        J_drz_dq3=diff(r_z,q3);
        J_drz_dq4=diff(r_z,q4);

        J_velocity_vector=[J_drx_dq1 J_drx_dq2 J_drx_dq3 J_drx_dq4; J_dry_dq1 J_dry_dq2 J_dry_dq3 J_dry_dq4 ; J_drz_dq1 J_drz_dq2 J_drz_dq3 J_drz_dq4];
    else
    end
    
    if q_initial(2)>15
        J_drx_dq1=0;
        J_drx_dq2=0;
        J_drx_dq3=diff(r_x,q3);
        J_drx_dq4=diff(r_x,q4);

        J_dry_dq1=0;
        J_dry_dq2=0;
        J_dry_dq3=diff(r_y,q3);
        J_dry_dq4=diff(r_y,q4);

        J_drz_dq1=0;
        J_drz_dq2=0;
        J_drz_dq3=diff(r_z,q3);
        J_drz_dq4=diff(r_z,q4);

        J_velocity_vector=[J_drx_dq1 J_drx_dq2 J_drx_dq3 J_drx_dq4; J_dry_dq1 J_dry_dq2 J_dry_dq3 J_dry_dq4 ; J_drz_dq1 J_drz_dq2 J_drz_dq3 J_drz_dq4];

    else
    end
    error_2=abs((P_2_total_norm-P_interp_2_total_norm)/P_2_total_norm)
    count=count+1;
% Uncomment to Plot Position 
    x_2(count,1)=P_interp_2_total(1,1);
    y_2(count,1)=P_interp_2_total(2,1);
    z_2(count,1)=P_interp_2_total(3,1);
   

% Uncomment to Plot Angular Velocity
%      theta_x_2(count,1)=P_interp_2_total(4,1);
%      theta_y_2(count,1)=P_interp_2_total(5,1);
%      theta_z_2(count,1)=P_interp_2_total(6,1);
%    

 end
 error_2
 
% Uncomment to Plot Position 
 plot3(x_2,y_2,z_2,'g-o')
title('Plotted trajectory of Position for Iteration 2');

% Uncomment to Plot Angular Velocity
% plot3(theta_x_2,theta_y_2,theta_z_2,'g-o')
% title('Plotted trajectory of Angle for Iteration 2');
% 
%  
