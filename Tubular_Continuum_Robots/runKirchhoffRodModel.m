clc
clear vars
close all
global alpha1
global alpha2
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
global g_global
global s_global;

% This script generates a struct tcr (Tubular Continuum Robot)
% which contains specific parameters of the robot, it executes the
% kirchhoff rod model and includes a function which can be used to visualize the robot.
%
%% Tube specifications
%Tube 1 (inner) parameters
L1s=(100)*.001; %straight length [mm]
L1c=(100)*.001; %curved length [mm]
ro1=1.75/2*.001;%m
ri1=1.25/2*.001;%m
u1s=[50;0;0]; %u1s constant initial curvature vector for the curved sections uis=[uisx;uisy;uisz]

I1=1/4*pi*(ro1^4-ri1^4);
J1=2*I1;
E1=60*10^9;%Pa
nu1=0.3;

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

[tcr] = setupTCR(L1s,L1c,u1s,ro1,ri1,E1,nu1,L2s,L2c,u2s,ro2,ri2,E1,nu1);

%configurational inputs
beta1=60*.001; %translation [mm]
beta2=20*.001;

alpha1=3.1416; %rotation 
alpha2=0;


%% this is where your kirchhoff model goes
u_ini=[0;0;0;0];
  options1=optimset('TolFun',1*exp(-6));
u=fsolve(@COMPUTE_RESIDUAL_ERROR,u_ini,options1);
g=g_global;
s=s_global;


%% visualize the robot
drawTCR(g,s,tcr,[-beta1,-beta2])

