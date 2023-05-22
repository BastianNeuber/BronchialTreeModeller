function M = FNC_CVCase48(Alpha,B,R_c)
% [Alpha(i+1);Alpha(i+2)],[B(j,:,1,i);B(j,:,2,i)],R_C(j,i)
%Die Funktion berechnet zu zwei gegebenen, sich schneidenden Halbgeraden den
%Mittelpunkt des Kreises, der beide Halbgeraden berührt und dessen z-Koordinate
%größer ist. Außerdem werden die beiden Berührpunkte mit den Halbgeraden
%ausgegeben.
%   Eingabewerte:
%       Alpha    Vektor der Winkel, den die Halbgeraden mit der z-Achse bilden.
%                Dabei ist immer der betragsmäßig kleinste
%                vorzeichenbehaftete (umgekehrter math. Sinn) gemeint.
%       B        Stürzpunkte der Geraden
%       Rc       Kreisradius
%       F        Fall

M = zeros(1,3);
epsilon1 = tan(pi/2-Alpha(1));
epsilon2 = B(1,3) - tan(pi/2-Alpha(1))*B(1,1) + R_c/sin(Alpha(1));
epsilon3 = tan(pi/2-Alpha(2));
epsilon4 = B(2,3) - tan(pi/2-Alpha(2))*B(2,1) - R_c/sin(Alpha(2));
M(1,1) = (epsilon4-epsilon2)/(epsilon1-epsilon3);
M(1,3) = epsilon1*M(1) + epsilon2;

end

