function [Pi,Pa,Pd] = FNC_Erzeugung_ESbd3D(Pi,Pa,Pd)

% Umformatieren
Pi = [Pi(:,1),zeros(size(Pi,1),1),Pi(:,2)];
Pa = [Pa(:,1),zeros(size(Pa,1),1),Pa(:,2)];

% Auflösung bestimmen
n = norm(Pi(1)-Pi(2));

% Innenfläche
Pi3d = [];
for i = 1 : size(Pi,1)
    N = ceil(2*pi*Pi(i,1)/n);
    for j = 1 : N-1
        Pi3d(end+1,:) = [cos(2*pi*j/N)*Pi(i,1),sin(2*pi*j/N)*Pi(i,1),Pi(i,3)];
    end
end
Pi = [Pi;Pi3d];

% Außenfläche
Pa3d = [];
for i = 1 : size(Pa,1)
    N = ceil(2*pi*Pa(i,1)/n);
    for j = 1 : N-1
        Pa3d(size(Pa3d,1)+1,:) = [cos(2*pi*j/N)*Pa(i,1),sin(2*pi*j/N)*Pa(i,1),Pa(i,3)];
    end
end
Pa = [Pa;Pa3d];

% Deckel
if isempty(Pd) == false
    Pd = [Pd(:,1),zeros(size(Pd,1),1),Pd(:,2)];
    Pd3d = [];
    for i = 1 : size(Pd,1)
        N = ceil(2*pi*Pd(i,1)/n);
        for j = 1 : N-1
            Pd3d(end+1,:) = [cos(2*pi*j/N)*Pd(i,1),sin(2*pi*j/N)*Pd(i,1),Pd(i,3)];
        end
    end
    Pd = [Pd;Pd3d];
end
end

