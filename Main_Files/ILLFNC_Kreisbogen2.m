function [P] = FNC_Kreisbogen2(M,Radius_oder_Punkt,Startwinkel,Winkel)
%Die Funktion gibt 100 Punkte aus, die auf einem Kreisbogen um den
%Mittelpunkt M mit dem Radius R von einem Startwinkel über einen Winkel
%liegen. Alternativ können auch ein Mittelpunkt, ein Punkt und ein Winkel
%angegeben werden.

%Falls Radius_oder_Punkt als Radius angegeben wurde
if size(Radius_oder_Punkt,1) == 1 && size(Radius_oder_Punkt,2) == 1
    R = Radius_oder_Punkt;
    alpha = linspace(Startwinkel,Startwinkel + Winkel);
    P = zeros(100,2);
    for i = 1 : 100
        P(i,:) = M + R*[cos(alpha(i)),sin(alpha(i))];
    end
%Falls Radius_oder_Punkt als Punkt angegeben wurde
else
    R = norm(M-Radius_oder_Punkt);
    Punkt = Radius_oder_Punkt;
    if M(1) < Punkt(1)
        Startwinkel = asin((Punkt(2)-M(2))/R);
    elseif M(1) > Punkt(1) && Punkt(2) < M(2)
        Startwinkel = -pi - asin((Punkt(2)-M(2))/R);
    elseif M(1) > Punkt(1) && Punkt(2) > M(2)
        Startwinkel = pi - asin((Punkt(2)-M(2))/R);
    end
    alpha = linspace(Startwinkel,Startwinkel + Winkel);
    P = zeros(100,2);
    for i = 1 : 100
        P(i,:) = M + R*[cos(alpha(i)),sin(alpha(i))];
    end
end
end

