function [K,x,y] = FNC_Winkelbogen(M,r,Startwinkel,Winkel,Pos,Xlim,Ylim)

K = FNC_Kreisbogen2(M,r,Startwinkel,Winkel);
[x,y] = FNC_AnnotCoord(Pos,Xlim,Ylim,K(99,:),K(100,:));

end