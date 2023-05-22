function [P] = FNC_Drehung3xz(A,B,alpha,P)
%Drehung3xz dreht eine Punktewolke mit dem Winkel alpha um eine beliebige
%Achse in der x-z-Ebene. Die Achse wird durch zwei Punkte A und B
%definiert.

% Verschiebung zum Ursprung
P = P - ones(size(P,1),1).*A;

% Drehung
P = FNC_Drehung3xz0(P,A-B,alpha);

%Rückverschiebung
P = P + ones(size(P,1),1).*A;

end

