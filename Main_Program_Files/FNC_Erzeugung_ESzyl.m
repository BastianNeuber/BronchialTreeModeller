function [P] = FNC_Erzeugung_ESzyl(n_r,L,Di,Da,Param_schliessen)

Pi = [];
Pd = [];
Pa = [];

%Innenzylinder
if L > 0  % Falls Abschluss-Rohrelemente erzeugt werden sollen...
    Pi = FNC_Erzeugung_RE(n_r,L,Di);
end

if Di ~= Da  % Falls eine Wandstärke berechnet wird...
    % Außenzylinder
    if L > 0  % Falls Abschluss-Rohrelemente erzeugt werden sollen...
        Pa = FNC_Erzeugung_RE(n_r,L,Da);
    end
    % Deckel
    if Param_schliessen == 1
        for i = ceil((Di/Da)*n_r) : n_r
            for j = 1 : n_r
                Pd = [Pd;[(Da/2)*(i/n_r)*cos(j/n_r*2*pi),(Da/2)*(i/n_r)*sin(j/n_r*2*pi),L]];
            end
        end
    end
elseif Param_schliessen == 1
    for j = 1 : n_r
        Pd = FNC_Erzeugung_Deckel(n_r,0,Da);
        Pd = Pd + [zeros(size(Pd,1),2),ones(size(Pd,1),1)*L];
    end
end

P = [Pi;Pa;Pd];

end

