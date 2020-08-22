function diskPoints = GeometricModel(TACR,q)
if abs(sum(q(1,:)))+abs(sum(q(2,:)))>1E-10
    msgbox('scheisse Daten')
else    
diskPoints=zeros(20,12);
rotation1=atan2(q(1,1)*cos(2*pi/3)-q(1,3),q(1,1)*sin(2*pi/3));
projection11=TACR.diskPitchRadius(1,1)*cos(rotation1);
theta1=pi/2+q(1,1)/projection11';
k1=((pi/2)-theta1)/TACR.segmentLength(1,1);
rotation2=atan2(q(2,1)*cos(2*pi/3)-q(2,3),q(2,1)*sin(2*pi/3));
projection21=TACR.diskPitchRadius(2,1)*cos(rotation2);
theta2=pi/2+q(2,1)/projection21';
k2=((pi/2)-theta2)/TACR.segmentLength(2,1);
s1=0:TACR.segmentLength(1,1)/(TACR.ndisks(1)-1):TACR.segmentLength(1,1);
i=2;
diskPoints(1,:)=[0,0,0,TACR.diskPitchRadius(1),0,0,-TACR.diskPitchRadius(1)/2,sqrt(3)*TACR.diskPitchRadius(1)/2,0,-TACR.diskPitchRadius(1)/2,-sqrt(3)*TACR.diskPitchRadius(1)/2,0];
while i<=10
 if k1==0
     Ti=[1,0,0,0;0,1,0,0;0,0,1,(i-1)*TACR.segmentLength(1,1)/9;0,0,0,1];
 else
  Ti=[cos(rotation1)*cos(k1*s1(i)),-sin(rotation1),cos(rotation1)*sin(k1*s1(i)),cos(rotation1)*(1-cos(k1*s1(i)))/k1;
    sin(rotation1)*cos(k1*s1(i)),  cos(rotation1), sin(rotation1)*sin(k1*s1(i)),sin(rotation1)*(1-cos(k1*s1(i)))/k1;
    -sin(k1*s1(i)),                0,              cos(k1*s1(i)),               sin(k1*s1(i))/k1;
    0,                             0,              0,                           1]*[cos(-rotation1),-sin(-rotation1),0,0;
       sin(-rotation1),cos(-rotation1), 0,0;
       0,              0,               1,0;
       0,              0,               0,1];
    
 end

A=Ti*[diskPoints(1,1:3),1]';
B=Ti*[diskPoints(1,4:6),1]';
C=Ti*[diskPoints(1,7:9),1]';
D=Ti*[diskPoints(1,10:12),1]';
diskPoints(i,1:3)=A(1:3);
diskPoints(i,4:6)=B(1:3);
diskPoints(i,7:9)=C(1:3);
diskPoints(i,10:12)=D(1:3);
 i=i+1;
end

s2=0:TACR.segmentLength(2,1)/TACR.ndisks(2,1):TACR.segmentLength(2,1);
i=11;
while i<=20
 if k2==0
      
     Tj=Ti*[1,0,0,0;0,1,0,0;0,0,1,(i-10)*TACR.segmentLength(2,1)/10;0,0,0,1];
      
   
 else
     Tj=Ti*[cos(rotation2)*cos(k2*s2(i-9)),-sin(rotation2),cos(rotation2)*sin(k2*s2(i-9)),cos(rotation2)*(1-cos(k2*s2(i-9)))/k2;
            sin(rotation2)*cos(k2*s2(i-9)),cos(rotation2), sin(rotation2)*sin(k2*s2(i-9)),sin(rotation2)*(1-cos(k2*s2(i-9)))/k2;
            -sin(k2*s2(i-9)),              0,              cos(k2*s2(i-9)),               sin(k2*s2(i-9))/k2;
            0,                             0,              0,                             1]*[cos(-rotation2),-sin(-rotation2),0,0;
       sin(-rotation2),cos(-rotation2), 0,0;
       0,              0,               1,0;
       0,              0,               0,1];
        
 end
 
v=[TACR.diskPitchRadius(2)/TACR.diskPitchRadius(1),TACR.diskPitchRadius(2)/TACR.diskPitchRadius(1),1,1];
V=diag(v);
A=Tj*(V*[diskPoints(1,1:3),1]');
B=Tj*(V*[diskPoints(1,4:6),1]');
C=Tj*(V*[diskPoints(1,7:9),1]');
D=Tj*(V*[diskPoints(1,10:12),1]');
diskPoints(i,1:3)=A(1:3);
diskPoints(i,4:6)=B(1:3);
diskPoints(i,7:9)=C(1:3);
diskPoints(i,10:12)=D(1:3);
 i=i+1;
end
end 


end