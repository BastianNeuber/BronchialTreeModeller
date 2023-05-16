function [R] = IFNC_Carinalradius(Param,Geom,j,n)
%IFNC berechnet einen Radius für eine Carinalverrundung auf der Basis der
%Geometrie von Tochterelementen.
%   Eingabewerte:
%       Param  Zeilenvektor mit Parametern
%       Geom   Matrix (2x2) mit Informationen über die Geometrie der
%              Tochterelemente: In der ersten Zeile stehen die Durchmesser
%              der Tochterelemente, in der zweiten Zeile stehen die
%              Verzweigungswinkel der Tochterelemente vom Mutterelement
%       j      Nummer der Punktreihe
%       n      Gesamtzahl der Punktreihen auf einem Rohrelement

R = sin(pi*(2*j/n))*(Param(1)*(Geom(1,1)+Geom(1,2)) + Param(2)*(Geom(2,1)+Geom(2,2)));

end

