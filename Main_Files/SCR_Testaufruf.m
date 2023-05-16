%% FNC_HeistracherBifurkation

clc
clear
hold off
close all

D = [8;4.5;2.5];
Alpha = [0;(50/180)*(pi);(-70/180)*(pi)];
sigma_M = [1,1.5,0,0,1]; % Konstruktionsmittelpunkte
sigma_S = [3,1.3]; % Sigmoidfunktion
n = 102;
sigma_CR = [0.1,0,0]; % Carinalverrundung
[P1] = FNC_Erzeugung_VE(D,Alpha,sigma_M,n,sigma_S(2),sigma_S(1),sigma_CR);
% sigma_CR = [0.0000001,0,0.5];
% [P2] = FNC_Erzeugung_VE(D,Alpha,sigma_M,n,sigma_S(2),sigma_S(1),sigma_CR);
% sigma_CR = [0.0000001,0,-0.3];
% [P3] = FNC_Erzeugung_VE(D,Alpha,sigma_M,n,sigma_S(2),sigma_S(1),sigma_CR);

% figure
scatter3(P1(:,1,1),P1(:,2,1),P1(:,3,1),2);
hold on
% scatter3(P2(:,1,1),P2(:,2,1),P2(:,3,1),2);
% scatter3(P3(:,1,1),P3(:,2,1),P3(:,3,1),2);
% xlim([-10,10]);
% ylim([-10,10]);
% zlim([0,12]);
view([0,0,-1]);

axis equal

%%
% 8°26'44.02'' E / 49°35'36.91'' N
E = [8,26,44.02];   % 8°57'42.47'' E
N = [49,35,36.91];   % 48°29'32.18'' N
Ed = E(1) + E(2)/60 + E(3)/3600;
Nd = N(1) + N(2)/60 + N(3)/3600;
[Ed, Nd]

%% Kollisionsprüfung

A = [4,3;5,2;6,1;7,0];
B = [1,2;4,1;0,3];

Delta = zeros(size(B,1),1);
for i = 1 : size(B,1)
    delta = (ones(size(A,1),size(B,2)).*B(i,:) - A);
    Delta(i) = min(sqrt(delta(:,1).^2+delta(:,2).^2));
end
min(Delta)

% Delta = zeros(size(A,1),1);
% for i = 1 : size(A,1)
%     Delta(i) = norm(A(i,:)-x);
% end
% Minimum = min(Delta)

%% Sockel

D = 25;
D0 = 4;
L = 50;
n = 50;
P = FNC_Erzeugung_Sockel(D0,L,D,n);
scatter3(P(:,1),P(:,2),P(:,3),2);
axis equal


%% Carinalauflösung

load('test.mat');


%% Ellipse
close all
a = 2;
b = 1;
c = 10;
d = 2;
x = linspace(-a+c,a+c);
yp = zeros(100,1);
ym = zeros(100,1);
for i = 1 : 100
    yp(i) = b/a*sqrt(a^2-(c-x(i))^2)+d;
end
for i = 1 : 100
    ym(i) = -b/a*sqrt(a^2-(c-x(i))^2)+d;
end
eta = linspace(0,pi);
for i = 1 : 100
    ye = sqrt(1-((sin(eta(i))*b)^2)/(b^2))*a;
end
plot(x(:),ye(:),'r');
% plot(x(:),yp(:),'b');
% hold on
% plot(x(:),ym(:),'b');
axis equal


%% Surf

[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
C = X.*Y;
surf(X,Y,Z,C)
colorbar

axis equal


%%
% z wurde in Abh. von x und y gemessen
% Festlegen des Gitters in x und y Koordinate
rangeX  = -10:0.1:10;
rangeY  = -10:0.1:10;
[X,Y]=meshgrid(rangeX,rangeY);
% Interpolation der Messwerte Z an den Gitterpunkten X,Y
Z=griddata(x,y,z,X,Y,'cubic');
% Plot als Fläche
surf(X,Y,Z)
hold off
grid on 




