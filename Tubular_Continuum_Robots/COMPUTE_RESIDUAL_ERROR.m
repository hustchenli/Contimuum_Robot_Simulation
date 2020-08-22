function error = COMPUTE_RESIDUAL_ERROR(u_ini)
global g_global;
global s_global;
global alpha1
global alpha2
global L1s
global L1c
global L2s
global L2c
global u1s
global I1
global J1
global E1
global u2s
global I2
global J2
global E2
global beta1
global beta2
global G1
global G2


u_ini=[0;0;0;0];
%% Tube specifications
%Tube 1 (inner) parameters

r_ini=[0;0;0];
w=alpha1;
R_ini=[cos(w),-sin(w),0    ;
       sin(w),cos(w) ,0    ;
       0    , 0      ,1     ];
theta1=0;
theta2=alpha2-alpha1;
s_ini=0;
y_ini=[(r_ini)',reshape(R_ini,1,9),(u_ini)',theta2,s_ini];
s=s_ini;
y=y_ini;
%pick up the positive transition points
T=[0,L1s-beta1,L2s-beta2,L1s+L1c-beta1,L2s+L2c-beta2];
Q=[-1,sort(T)];
p=[];
pen=1;
for dic=2:length(T)+1
    if Q(dic)>=0&&Q(dic)>Q(dic-1)
        p(pen)=Q(dic);
        pen=pen+1;
    end
end

for j=1:length(p)-1 
  options2=odeset('initialStep',.001);
  [a,b]= ode45(@(s,y) COMPUTE_DERIVATIVE(s,y),[p(j),p(j+1)],y_ini,options2);
  s=[s;a];
  y=[y;b];
 
 
  
  if j~=length(p)-1
     
      u1_bef=[y(end,13);y(end,14);y(end,15)];
      R1=[cos(alpha1),-sin(alpha1),0;sin(alpha1),cos(alpha1),0;0,0,1];
      R_theta1=[cos(theta1) -sin(theta1) 0
                sin(theta1) cos(theta1)  0
                0           0            1];
      R_theta2=[cos(theta2) -sin(theta2) 0
                sin(theta2) cos(theta2)  0
                0           0            1];
      R2=R1*R_theta2;
      
     
      K1=[E1*I1,0,0;0,E1*I1,0;0,0,G1*J1];
      K2=[E2*I2,0,0;0,E2*I2,0;0,0,G2*J2];
      
     %determin u and k before and after each transition point 
    if (s(end)<L1s-beta1)                                                          
    u1_star_bef=[0;0;0];
    u1_star_aft=[0;0;0];
    K1_bef=K1;
    K1_aft=K1;
    elseif (s(end)>L1s-beta1&&s(end)<L1s+L1c-beta1)                   
    u1_star_bef=u1s;
    u1_star_aft=u1s;
    K1_bef=K1;
    K1_aft=K1;
    elseif (s(end)==L1s-beta1)                           
    u1_star_bef=[0;0;0];                         
    u1_star_aft=u1s; 
    K1_bef=K1;
    K1_aft=K1;
    elseif (s(end)==L1s+L1c-beta1)
    u1_star_bef=u1s;                         
    u1_star_aft=[0;0;0];
    K1_bef=K1;
    K1_aft=zeros(3,3);
    elseif (s(end)>L1s+L1c-beta1)
         u1_star_bef=[0;0;0];
         u1_star_aft=[0;0;0]; 
         K1_bef=zeros(3,3);
         K1_aft=zeros(3,3);
    end
if (s(end)<L2s-beta2)                               
    u2_star_bef=[0;0;0];
    u2_star_aft=[0;0;0];
    K2_bef=K2;
    K2_aft=K2;
elseif (s(end)>L2s-beta2&&s(end)<L2s+L2c-beta2)                   
    u2_star_bef=u2s;
    u2_star_aft=u2s;
    K2_bef=K2;
    K2_aft=K2;
elseif (s(end)==L2s-beta2)                          
    u2_star_bef=[0;0;0];                        
    u2_star_aft=u2s;
    K2_bef=K2;
    K2_aft=K2;
elseif (s(end)==L2s+L2c-beta2)
    u2_star_bef=u2s;                        
    u2_star_aft=[0;0;0];
    K2_bef=K2;
    K2_aft=zeros(3,3);
elseif (s(end)>L2s+L2c-beta2) 
    u2_star_bef=[0;0;0];
    u2_star_aft=[0;0;0];
    K2_bef=zeros(3,3);
    K2_aft=zeros(3,3);
end



     
     u1z_bef=y(end,15);
      u2z_bef=y(end,16);
     u1z_aft=u1_star_aft(3)+u1z_bef-u1_star_bef(3);
     u2z_aft=u2_star_aft(3)+u2z_bef-u2_star_bef(3);
    
      R1=reshape(y(end,4:12),3,3);
      R2=R1*R_theta2;
  
   u1_aft=(((R1*K1_aft+R2*K2_aft*R_theta2')))\(R1*(K1_bef*u1_bef-K1_bef*u1_star_bef+K1_aft*u1_star_aft)+R2*(K2_bef*R_theta2'*u1_bef-K2_bef*u2_star_bef+K2_aft*u2_star_aft));
      u1x_aft=u1_aft(1);
      u1y_aft=u1_aft(2);
      
     
     
      
      y_ini=[y(end,1:12),u1x_aft,u1y_aft,u1z_aft,u2z_aft,y(end,17:end)];

  end
  
end
%Apply the boundary conditions to compute the residual error.
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
error=y(end,13:16)'-[u1_star;u2_star(3)];
[z,s_palt]=size(y);
c=ones(z,1);
g_global=zeros(z,16);
g_global(:,1:3)=y(:,4:6) ;
g_global(:,5:7)=y(:,7:9);
g_global(:,9:11)=y(:,10:12);
g_global(:,13:15)=y(:,1:3);
g_global(:,16)=c;
s_global=s;
end
      
  
