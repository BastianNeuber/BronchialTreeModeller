function [P,M] = FNC_Box(Pos,B,H)

Pos = [Pos(1),Pos(2)];
P = FNC_Punkte2Strecke([-B/2,H/2],[B/2,H/2],99);
M(1,:) = P(49,:);
P = [P;FNC_Punkte2Strecke([B/2,H/2],[B/2,-H/2],99)];
M(2,:) = P(148,:);
P = [P;FNC_Punkte2Strecke([B/2,-H/2],[-B/2,-H/2],99)];
M(2,:) = P(247,:);
P = [P;FNC_Punkte2Strecke([-B/2,-H/2],[-B/2,H/2],99)];
M(2,:) = P(346,:);
P = P + ones(size(P,1),1).*Pos;
M = M + ones(size(M,1),1).*Pos;

end

