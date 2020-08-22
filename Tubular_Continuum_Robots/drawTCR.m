function [fig] = drawTCR(g,s,tcr,b)
%% drawTCR.m
% This function visualizes a tubular continuum robot with n tubes
%
%
% input: g(i,:) = [R(1,1) R(1,2) R(1,3) 0 R(2,1) R(2,2) R(2,3) 0 R(3,1) R(3,2) R(3,3) 0 r(1) r(2) r(3) 1];
%        with rotation matrix R and position vector r
%        s is the arc length variable
%        struct tcr
%        b: beta1 and beta2, translational inputs 

    %transparency factor
    alpha=1; %1 = opaque, 0 = completely transparent
    fig = gcf;
        
    numberTubes = tcr;   
    color = zeros(length(numberTubes),3);
    radius = zeros(length(numberTubes),1);
    
    k = 0.4; 
    for i = 1:length(numberTubes)
        %color of tubes
        color(i,:) = [k k k];
        k = k + 0.2;
        %radius
        radius(i) = tcr(i).tube.ro;
    end
    
    %determine tube_ends
    r = ones(length(s),1)*tcr(1).tube.ro;
    for i=1:length(s)            
        for j=1:length(numberTubes)
             if (s(i) >= b(j)) && (s(i) <= b(j) + tcr(j).tube.L)
                        r(i)=tcr(j).tube.ro;
             end    
        end            
    end
        
    plotsection(g,radius(1),color(1,:),alpha);
        
    for j=2:length(numberTubes)
        indeces=find(r==tcr(j).tube.ro);
        plotsection(g(indeces(1):indeces(length(indeces)),:),radius(j),color(j,:),alpha);
    end
       

    hold on
    color_entry_point = [1 1 1]*0.9;
    squaresize = 0.02;
    thickness = 0.001;
    patch([-1 1 1 -1]*squaresize,[-1 -1 1 1]*squaresize,[-1 -1 -1 -1]*thickness,color_entry_point)
    patch([-1 1 1 -1]*squaresize,[-1 -1 1 1]*squaresize,[0 0 0 0],color_entry_point)
    patch([1 1 1 1]*squaresize,[-1 -1 1 1]*squaresize,[-1 0 0 -1]*thickness,color_entry_point)
    patch([-1 -1 -1 -1]*squaresize,[-1 -1 1 1]*squaresize,[-1 0 0 -1]*thickness,color_entry_point)
    patch([-1 1 1 -1]*squaresize,[-1 -1 -1 -1]*squaresize,[-1 -1 0 0]*thickness,color_entry_point)
    patch([-1 1 1 -1]*squaresize,[1 1 1 1]*squaresize,[-1 -1 0 0]*thickness,color_entry_point)

    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    grid on
    view([0.5 0.5 0.5])
    daspect([1 1 1])
end
