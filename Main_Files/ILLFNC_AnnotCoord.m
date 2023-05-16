function [x,y] = FNC_AnnotCoord(Pos,Xlim,Ylim,A,B)
%Die Funktion berechnet aus get(gca,'position'), get(gca,'xlim'),
%get(gca,'ylim') und zwei Koordinatenpunkten die beiden Parameter x und y,
%die für das zeichnen eines Pfeils mithilfe der annotation-Funktion
%erstellt werden können. mode gibt an, ob ein Doppel- oder Einfachpfeil
%('d' oder 's') gezeichnet wird.

x = zeros(1,2);
y = x;

x(1) = Pos(1)+(A(1)-Xlim(1))/(Xlim(2)-Xlim(1)) * Pos(3);
x(2) = Pos(1)+(B(1)-Xlim(1))/(Xlim(2)-Xlim(1)) * Pos(3);
y(1) = Pos(2)+(A(2)-Ylim(1))/(Ylim(2)-Ylim(1)) * Pos(4);
y(2) = Pos(2)+(B(2)-Ylim(1))/(Ylim(2)-Ylim(1)) * Pos(4);

end

