
function dy = COMPUTE_DERIVATIVE(s,y)


global L1s
global L1c
global L2s
global L2c
global ro1
global ri1
global u1s
global I1
global J1
global E1
global nu1
global ro2
global ri2
global u2s
global I2
global J2
global E2
global nu2
global beta1
global beta2
global G1
global G2

e3=[0;0;1];
%theta1=0;
theta2=y(17);
%Tube 1 (inner) parameters
L1s=(100)*.001; %straight length [mm]
L1c=(100)*.001; %curved length [mm]
ro1=1.75/2*.001;%m
ri1=1.25/2*.001;%m
u1s=[50;0;0]; %u1s constant initial curvature vector for the curved

I1=1/4*pi*(ro1^4-ri1^4);
J1=2*I1;
E1=60*10^9;%Pa
nu1=0.3;
G1=E1/(2*(1+nu1));

%Tube 2 (outer) parameters
L2s=(50)*.001;
L2c=(50)*.001;
ro2=2.50/2*.001;%m
ri2=2.00/2*.001;%m
u2s=[10;0;0];

I2=1/4*pi*(ro2^4-ri2^4);
J2=2*I2;
E2=60*10^9;%Pa
nu2=0.3;
G2=E2/(2*(1+nu2));
%Stiffness matrix
K1=[E1*I1,0,0;0,E1*I1,0;0,0,G1*J1];
K2=[E2*I2,0,0;0,E2*I2,0;0,0,G2*J2];

%determin u and k before and after each transition point 
if s(end) >=(L2s+L2c-beta2)
    K2=zeros(3,3);
end

K=K1+K2;
if s(end)>(L1s+L1c-beta1)
    K=zeros(3,3);
end
if s(end)<=L1s-beta1
          u1_star=[0;0;0];
      elseif (s(end)>L1s-beta1&&s(end)<=L1s+L1c-beta1)
          u1_star=u1s;
      elseif (s(end)>L1s+L1c-beta1)
         
          u1_star=[0;0;0];
 end
   if s(end)<=L2s-beta2
          u2_star=[0;0;0];
      elseif (s(end)>L2s-beta2&&s(end)<=L2s+L2c-beta2)
          u2_star=u2s;
      elseif (s(end)>L2s+L2c-beta2)
         
          u2_star=[0;0;0];
    end
  
 
%Differential Equations

R_theta2=[cos(theta2) -sin(theta2) 0
          sin(theta2) cos(theta2)  0
          0           0            1];
u1x=y(13);
u1y=y(14);
u1z=y(15);
u2z=y(16);
u1=[u1x;u1y;u1z];
theta2_dot=u2z-u1z;
u2=(R_theta2)'*u1+theta2_dot*[0;0;1];
u2x=u2(1);
u2y=u2(2);
u2z=u2(3);
R_theta2_dot=[-sin(theta2) -cos(theta2) 0
          cos(theta2) -sin(theta2)  0
          0           0            1];

u1_dach=[0,-u1z,u1y;u1z,0,-u1x;-u1y,u1x,0];
u2_dach=[0,-u2z,u2y;u2z,0,-u2x;-u2y,u2x,0];
u1_star_dot=[0;0;0];

u2_star_dot=[0;0;0];
u1xy_dot=-inv(K)*(K1*(-u1_star_dot)+u1_dach*K1*(u1-u1_star)+R_theta2*(K2*(theta2_dot*R_theta2_dot*u1-u2_star_dot)+u2_dach*K2*(u2-u2_star)));

u1x_dot=u1xy_dot(1);
u1y_dot=u1xy_dot(2);
u1z_dot=u1_star_dot(3)+(E1*I1)/(G1*J1)*(u1x*u1_star(2)-u1y*u1_star(1));
u2z_dot=u2_star_dot(3)+(E2*I2)/(G2*J2)*(u2x*u2_star(2)-u2y*u2_star(1));
R_dot=reshape(y(4:12),3,3)*u1_dach;
r_dot=reshape(y(4:12),3,3)*e3;
s_dot=1;
dy=[r_dot;reshape(R_dot,[9,1]);u1x_dot;u1y_dot;u1z_dot;u2z_dot;theta2_dot;s_dot];

end