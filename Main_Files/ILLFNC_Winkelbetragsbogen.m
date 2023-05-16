function [K,x1,y1,x2,y2] = FNC_Winkelbetragsbogen(M,r,Startwinkel,Winkel,Pos,Xlim,Ylim)

K = FNC_Kreisbogen2(M,r,Startwinkel,Winkel);
[x1,y1] = FNC_AnnotCoord(Pos,Xlim,Ylim,K(99,:),K(100,:));
[x2,y2] = FNC_AnnotCoord(Pos,Xlim,Ylim,K(2,:),K(1,:));

end