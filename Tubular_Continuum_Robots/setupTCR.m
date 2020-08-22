function [tcr] = setupTCR(L1s,L1c,u1s,ro1,ri1,E1,nu1,L2s,L2c,u2s,ro2,ri2,E2,nu2)
%% drawTCR.m
% This function creates a stucture tcr
%
%
% Creates a struct ctcr 
% n=2 tubes, with tube 1 the innermost tube
% Lns: length of the straight section (m)
% Lnc: length of the curved section (m)
% u1sxn: curvature of the curved section (1/m)
% ron,rin: radius of the outer and inner tube (m)
% material properties E, nu 

    G1 = E1/(2*(1+nu1));
    G2 = E2/(2*(1+nu2));    
    
    %Tube 1, inner
    tube1_I = 1/4*pi*(ro1^4-ri1^4);
    tube1 = struct('L',L1s+L1c,'Ls',L1s,'Lc',L1c,...
                   'ro',ro1,'ri',ri1,...
                   'u1s',u1s,'r',1/u1s,...
                   'I',tube1_I,'E',E1,'G',G1,'J',2*tube1_I);
               
    %Tube 2, outer
    tube2_I = 1/4*pi*(ro2^4-ri2^4);
    tube2 = struct('L',L2s+L2c,'Ls',L2s,'Lc',L2c,...
                   'ro',ro2,'ri',ri2,...
                   'u2s',u2s,'r',1/u2s,...
                   'I',tube2_I,'E',E2,'G',G2,'J',2*tube2_I);

    
    tcr(1) = struct('tube',tube1);
    tcr(2) = struct('tube',tube2);
        
end