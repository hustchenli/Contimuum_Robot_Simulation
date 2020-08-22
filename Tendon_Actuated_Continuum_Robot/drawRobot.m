function [] = drawRobot(diskPoints,TACR)

% This function visualizes a 2 tendon actuated continuum robot
%
% input: struct TACR, 
%        diskPoints: n rows for n disks, 12 columns for coordinates
%        (x,y,z) for central backbone, points tendon 1 (x,y,z), points tendon 2 (x,y,z),
%        points tendon 3 (x,y,z)
 
figure();
hold on;  

color = [148/255 148/255 148/255];

ndisks = TACR.ndisks;                   
diskRadius = TACR.diskRadius;            
diskHeight = TACR.diskHeight;            
diskPitchRadius = TACR.diskPitchRadius;
radialPoints = 20;  % the number of axial rectangles to draw on the cannula surface
circle_bot = zeros(3,radialPoints);
circle_top = zeros(3,radialPoints);
k=1;

for i=1:sum(ndisks)
    if i==ndisks(1,1)+1
        k=2;
    end
    n = ((cross(-diskPoints(i,1:3)+diskPoints(i,7:9),-diskPoints(i,1:3)...
        +diskPoints(i,4:6)))')/norm((cross(-diskPoints(i,1:3)+...
        diskPoints(i,7:9),-diskPoints(i,1:3)+diskPoints(i,4:6)))');
    
    for j=1:radialPoints
        rot_angle = 2*pi - j*(2*pi/radialPoints);
        n1 = n(1,1);
        n2 = n(2,1);
        n3 = n(3,1);
        R_n = [n1^2*(1-cos(rot_angle))+cos(rot_angle), n1*n2*(1-cos(rot_angle))-n3*sin(rot_angle), n1*n3*(1-cos(rot_angle))+n2*sin(rot_angle);...
            n2*n1*(1-cos(rot_angle))+n3*sin(rot_angle), n2^2*(1-cos(rot_angle))+cos(rot_angle), n2*n3*(1-cos(rot_angle))-n1*sin(rot_angle);...
            n3*n1*(1-cos(rot_angle))-n2*sin(rot_angle), n3*n2*(1-cos(rot_angle))+n1*sin(rot_angle), n3^2*(1-cos(rot_angle))+cos(rot_angle)];        
        circle_bot(:,j) = diskPoints(i,1:3)' + (R_n * (((-diskPoints(i,1:3)+diskPoints(i,4:6))/norm((-diskPoints(i,1:3)+diskPoints(i,4:6)))) * diskRadius(k,1))');
        circle_top(:,j) = circle_bot(:,j) + n*diskHeight;
    end
    
    %draw circle
    fill3(circle_bot(1,:),circle_bot(2,:),circle_bot(3,:),color);
    fill3(circle_top(1,:),circle_top(2,:),circle_top(3,:),color);
    
    %draw cylinder shell
    x=[circle_bot(1,1:end); circle_top(1,1:end)];
    y=[circle_bot(2,1:end); circle_top(2,1:end)];
    z=[circle_bot(3,1:end); circle_top(3,1:end)];
    surf(x,y,z,'FaceColor',color,'MeshStyle','row');
    
    % plot tendons of the 2nd segment through the 1st segment
    if i >= 1 && i <= ndisks(1,1) 
        t_1 = [diskPoints(i,1:3)+(-diskPoints(i,1:3)+diskPoints(i,4:6))/norm(-diskPoints(i,1:3)+diskPoints(i,4:6))*diskPitchRadius(2,1),...
            diskPoints(i+1,1:3)+(-diskPoints(i+1,1:3)+diskPoints(i+1,4:6))/norm(-diskPoints(i+1,1:3)+diskPoints(i+1,4:6))*diskPitchRadius(2,1)];
        t_2 = [diskPoints(i,1:3)+(-diskPoints(i,1:3)+diskPoints(i,7:9))/norm(-diskPoints(i,1:3)+diskPoints(i,7:9))*diskPitchRadius(2,1),...
            diskPoints(i+1,1:3)+(-diskPoints(i+1,1:3)+diskPoints(i+1,7:9))/norm(-diskPoints(i+1,1:3)+diskPoints(i+1,7:9))*diskPitchRadius(2,1)];
        t_3 = [diskPoints(i,1:3)+(-diskPoints(i,1:3)+diskPoints(i,10:12))/norm(-diskPoints(i,1:3)+diskPoints(i,10:12))*diskPitchRadius(2,1),...
            diskPoints(i+1,1:3)+(-diskPoints(i+1,1:3)+diskPoints(i+1,10:12))/norm(-diskPoints(i+1,1:3)+diskPoints(i+1,10:12))*diskPitchRadius(2,1)];
        plot3([t_1(1,1) t_1(1,4)],[t_1(1,2) t_1(1,5)],[t_1(1,3) t_1(1,6)],'LineWidth',1.1,'Color',[0 0 0]);
        plot3([t_2(1,1) t_2(1,4)],[t_2(1,2) t_2(1,5)],[t_2(1,3) t_2(1,6)],'LineWidth',1.1,'Color',[0 0 0]);
        plot3([t_3(1,1) t_3(1,4)],[t_3(1,2) t_3(1,5)],[t_3(1,3) t_3(1,6)],'LineWidth',1.1,'Color',[0 0 0]);
    end
end

%% plot backbone + tendons %%
tb_x = [diskPoints(1:sum(ndisks)-1,1) diskPoints(2:sum(ndisks),1)];
tb_y = [diskPoints(1:sum(ndisks)-1,2) diskPoints(2:sum(ndisks),2)];
tb_z = [diskPoints(1:sum(ndisks)-1,3) diskPoints(2:sum(ndisks),3)];
plot3(tb_x,tb_y,tb_z,'LineWidth',3,'Color',[28/255 28/255 28/255]);
t11_x = [diskPoints(1:ndisks(1,1)-1,4) diskPoints(2:ndisks(1,1),4)];
t11_y = [diskPoints(1:ndisks(1,1)-1,5) diskPoints(2:ndisks(1,1),5)];
t11_z = [diskPoints(1:ndisks(1,1)-1,6) diskPoints(2:ndisks(1,1),6)];
plot3(t11_x,t11_y,t11_z,'LineWidth',1.1,'Color',[0 0 0]);
t12_x = [diskPoints(1:ndisks(1,1)-1,7) diskPoints(2:ndisks(1,1),7)];
t12_y = [diskPoints(1:ndisks(1,1)-1,8) diskPoints(2:ndisks(1,1),8)];
t12_z = [diskPoints(1:ndisks(1,1)-1,9) diskPoints(2:ndisks(1,1),9)];
plot3(t12_x,t12_y,t12_z,'LineWidth',1.1,'Color',[0 0 0]);
t13_x = [diskPoints(1:ndisks(1,1)-1,10) diskPoints(2:ndisks(1,1),10)];
t13_y = [diskPoints(1:ndisks(1,1)-1,11) diskPoints(2:ndisks(1,1),11)];
t13_z = [diskPoints(1:ndisks(1,1)-1,12) diskPoints(2:ndisks(1,1),12)];
plot3(t13_x,t13_y,t13_z,'LineWidth',1.1,'Color',[0 0 0]);
t21_x = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,4) diskPoints(ndisks(1,1)+2:sum(ndisks),4)];
t21_y = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,5) diskPoints(ndisks(1,1)+2:sum(ndisks),5)];
t21_z = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,6) diskPoints(ndisks(1,1)+2:sum(ndisks),6)];
plot3(t21_x,t21_y,t21_z,'LineWidth',1.1,'Color',[0 0 0]);
t22_x = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,7) diskPoints(ndisks(1,1)+2:sum(ndisks),7)];
t22_y = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,8) diskPoints(ndisks(1,1)+2:sum(ndisks),8)];
t22_z = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,9) diskPoints(ndisks(1,1)+2:sum(ndisks),9)];
plot3(t22_x,t22_y,t22_z,'LineWidth',1.1,'Color',[0 0 0]);
t23_x = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,10) diskPoints(ndisks(1,1)+2:sum(ndisks),10)];
t23_y = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,11) diskPoints(ndisks(1,1)+2:sum(ndisks),11)];
t23_z = [diskPoints(ndisks(1,1)+1:sum(ndisks)-1,12) diskPoints(ndisks(1,1)+2:sum(ndisks),12)];
plot3(t23_x,t23_y,t23_z,'LineWidth',1.1,'Color',[0 0 0]);

axis equal; axis tight;
xlabel('x-axis [mm]');
ylabel('y-axis [mm]');
zlabel('z-axis [mm]');
view([-20 14]);
% grid on;

end