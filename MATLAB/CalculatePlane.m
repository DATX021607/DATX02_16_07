function [zplane] = CalculatePlane( P1, P2, P3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

normal = cross(P1-P3,P1-P2)
b1 = normal(1);
b2 = normal(2);
b3 = normal(3);

%A = 0 +(b1*B)+b2*C;
%plot3(B,C,A)
hold on



point2 = [b1 b2 b3];

pointsToPlot = [point2;P2];

b1 = b1/-b3
b2 = b2/-b3

%TODO: Make these values more general.
Xn = linspace(0,300);
Yn = linspace(20,600);
Zn = b1*Xn + b2*Yn;
hold on
%plot3(pointsToPlot(:,1),pointsToPlot(:,2), pointsToPlot(:,3))
hold on

syms x y z
P= [x y z];
planefunction = dot(normal, P-P1)
hold on
zplane = solve(planefunction, z);

end

