function [Xaxis,Yaxis] = FNC_CoordSys2(O,Pos,xlim,ylim,Xaxislim,Yaxislim)
%FNC_CoordSys2 erstellt die Koordinaten für zwei Pfeile, die im Plot als
%Koordinatensystem verwendet weden können.

Xaxis = zeros(2,2);
[Xaxis(1,:),Xaxis(2,:)] = FNC_AnnotCoord(Pos,xlim,ylim,[Xaxislim(1),O(2)],[Xaxislim(2),O(2)]);
Yaxis = zeros(2,2);
[Yaxis(1,:),Yaxis(2,:)] = FNC_AnnotCoord(Pos,xlim,ylim,[O(1),Yaxislim(1)],[O(1),Yaxislim(2)]);

end

