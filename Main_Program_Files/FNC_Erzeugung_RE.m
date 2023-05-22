function [P] = FNC_Erzeugung_RE(n,L,D)

n_KQ = FNC_Resolution2(n*L/(D*pi));

% Vorbereitung
n_KQ = ceil(n_KQ);
P = zeros(n*n_KQ,3);

% Erzeugung des Rohrelementes
for k = 1 : n_KQ
    for p = 1 : n
        P(p+(k-1)*n,:) = [(D/2)*cos(2*pi*p/n),(D/2)*sin(2*pi*p/n),(k-1)*(L/(n_KQ-1))];
    end
end