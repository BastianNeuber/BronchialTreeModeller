function [P] = FNC_Anpassung_bdES(n,Di,Da,Pi,Pa,Pd)

% Di    Innendurchmesser
% Da    Außendurchmesser
% Pf    Punkte der Form (2D)

%% Vorprüfung
if isempty(Pd) % Falls kein Deckel konstruiert wurde...
    Pd = [0,0,0];
end
if max([max(Pd(:,1)),max(Pa(:,1)),max(Pi(:,1))]) < Da/2  %  isempty(Pd) == false && 
    P = [];
else
    %% Punkte aus der Form löschen

    % Innen
    nPi = size(Pi,1);
    del = [];
    i = 1;
    Pi = sortrows(Pi,3);
    while i < nPi && sqrt(Pi(i,1)^2+Pi(i,2)^2) < Di/2
        del = [del;i];
        i = i + 1;
    end
    Pi(del,:) = [];

    % Außen
    nPa = size(Pa,1);
    del = [];
    i = 1;
    Pa = sortrows(Pa,3);
    while i < nPa && sqrt(Pa(i,1)^2+Pa(i,2)^2) < Da/2
        del = [del;i];
        i = i + 1;
    end
    Pa(del,:) = [];

    % Zusammenfassen
    P = [Pi;Pa;Pd];


    %% Nachbereitung

    % Verschiebung in z-Richtung
    DeltaA = min(Pa(:,3));
    DeltaI = min(Pi(:,3));
    P = P-[zeros(size(P,1),2),ones(size(P,1),1)]*DeltaA;

    % Verlängerung der Innenfläche durch Zylinder
    P = [P;FNC_Erzeugung_RE(n,DeltaI-DeltaA,2*Pi(1,1))];
end

end

