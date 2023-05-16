close
clc
clear

load 'Drehung.mat'

ABx = linspace(A(1),B(1));
ABy = linspace(A(2),B(2));
ABz = linspace(A(3),B(3));

plot3(ABx(:),ABy(:),ABz(:));
hold on
scatter3(P_original(:,1),P_original(:,2),P_original(:,3),1);

axis equal




%%
% n = [B(2) - A(2); -B(1) + A(1); 0];
% if norm(n) ~= 0
%     n = n/norm(n);
% 
%     %% Rotationswinkel
%     if B(3) > A(3)
%         alpha = asin(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A));
%     elseif B(1) > A(1)
%         alpha = (pi/2 + acos(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A)))*sign(B(1)-A(1));
%     elseif B(1) < A(1)
%         alpha = (pi + asin(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A)))*sign(B(1)-A(1));
%     end
%     s = sin(alpha); 
%     c = cos(alpha); 
% 
%     %% 3D rotation matrix:
%     x  = n(1);
%     y  = n(2);
%     z  = n(3);
%     mc = 1 - c;
%     R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
%         x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
%         x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];
% 
%     %% Rotation
%     P = P*R;