function [P] = FNC_Erzeugung_Sockel(D0,L,D,n,B)
% FNC_Erzeugung_Sockel erzeugt einen zylindrischen Sockel und einen runden
% Übergang zwischen dem Bronchialbaum und dem Sockel. Das Modell kann so
% ohne weitere konstruktive Schritte direkt in eine .stl umgewandelt
% werden.

n_Z = FNC_Resolution2(n*L/(D*pi));
P = zeros(n_Z*n+2*n^2,3);

% Übergang
for i = 1 : 4 : n
    for j = 1 : n
       P((i-1)*n+j,:) = [cos(j/n*2*pi)*(D/2-cos(i/n*pi/2)*(D-D0)/2),sin(j/n*2*pi)*(D/2-cos(i/n*pi/2)*(D-D0)/2),-sin(i/n*pi/2)*(D-D0)/2];
    end
end

% Zylinder
P(n^2+1 : n^2 + n*n_Z,:) = FNC_Erzeugung_RE(n,L,D)-[zeros(n_Z*n,1),zeros(n_Z*n,1),ones(n_Z*n,1)*(L+(D-D0)/2)];

% Loch bei Berechnung mit Wandstärke
if B ~= 0
    P = [P;-FNC_Erzeugung_RE(FNC_Resolution2(B/D*n),L+(D-D0)/2,B)];
end

% Boden
for i = ceil((B/D)*n) : 4 : n
    for j = 1 : n
        P = [P;[D/2*(i/n)*cos(j/n*2*pi),D/2*(i/n)*sin(j/n*2*pi),-(L+(D-D0)/2)]];
    end
end

end