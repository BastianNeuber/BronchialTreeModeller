function [dist,koll] = FNC_Abstandsberechnung(mode,A,B,D_A,D_B,Param_D,Param_dist)
% A           Punkte, die optimiert werden
% B           Fixpunkte
% D_A         Durchmesser aller Bronchien zu den Punkten aus A
% D_B         Durchmesser aller Bronchien zu den Punkten aus B
% Param_D     Gewichtung der Durchmesser
% Param_dist  Gewichtung der Abstände

% Abstand
if mode == 1
    dist = 0;
elseif mode == 2
    dist = Inf;
end

koll = [];  % Kollisionen

for i = 1 : size(A,1)
    for j = 1 : size(B,1)
        % Abstandsberechnung
        if mode == 1
            dist = dist + norm(A(i,:)-B(j,:))^Param_dist * (D_A(i)+D_B(j))^Param_D;
        elseif mode == 2 && norm(A(i,:)-B(j,:)) < dist
            dist = norm(A(i,:)-B(j,:));
        end
%         % Kollisionsprüfung
%         if norm(A(i,:)-B(j,:)) <= (D_A(i) + D_B(j))/2
%             koll = [koll;i,j];
%         end
    end
end
Halt = 1;