function [Mx,Mz] = FNC_CVCase2367(Alpha,B,R_c,R_k,M_k,CaseNum)
%FNC_CVCase2367 berechnet die Schnittpunkte eines Kreises mit einer Geraden.
%Dabei wird der betragsmäßig kleinere Wert als Lösung ausgegeben.

if CaseNum == 2
    RcFaktor = [-1,1];
elseif CaseNum == 3
    RcFaktor = [1,1];
elseif CaseNum == 6
    if Alpha > 0
        RcFaktor = [-1,1];
    else
        RcFaktor = [-1,-1];
    end
elseif CaseNum == 7
    if Alpha > 0
        RcFaktor = [1,-1];
    else
        RcFaktor = [-1,1];
    end
end

%Steigung der Geraden
epsilon1 = tan(pi/2-Alpha);

%z-Verschiebung der Geraden
epsilon2 = B(3)-tan(pi/2-Alpha)*B(1)+RcFaktor(1)*R_c/sin(Alpha);

%p-Anteil des Polynoms (vgl. p-q-Formel)
epsilon3 = 2*(epsilon1*epsilon2-M_k)/(1+epsilon1^2);

%q-Anteil des Polynoms (vgl. p-q-Formel)
epsilon4 = (M_k^2+epsilon2^2-(R_k+RcFaktor(2)*R_c)^2)/(1+epsilon1^2);

%Einsetzen der Lösung in die Geradengleichung
Mx = -epsilon3/2-sign(M_k)*sqrt((epsilon3^2)/4-epsilon4);
Mz = epsilon1*Mx+epsilon2;

end