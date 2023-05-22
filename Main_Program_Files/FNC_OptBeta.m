function [OptBetaA,OptBetaB,DMdouble,M_A,M_B] = FNC_OptBeta(D,M,DMdouble,Param_KV)


%% Position der Mutterverzweigung in M
empty = cellfun('isempty',M);  % Definition einer Cellfunction zur Ermittlung leerer Zellen
Gen = size(M,1)-sum(empty(:,1)==1);  % aktuelle Generation
nVerzw = size(M,2)-sum(empty(Gen,:)==1);
Pos = [Gen-1,nVerzw/2;Gen,nVerzw-1;Gen,nVerzw];
% A = [1. Tochterast außen
%      1. Tochterast innen
%      2. Tochterast außen
%      2. Tochterast innen
%      Mutterast unten
%      Mutterast oben
% B analog
A = [M{Pos(2,1),Pos(2,2)}{2,1}(1,:),D{Pos(2,1),Pos(2,2)}(1);...
     M{Pos(2,1),Pos(2,2)}{1,1}(1,:),D{Pos(2,1),Pos(2,2)}(1);...
	 M{Pos(2,1),Pos(2,2)}{2,1}(2,:),D{Pos(2,1),Pos(2,2)}(2);...
	 M{Pos(2,1),Pos(2,2)}{1,1}(2,:),D{Pos(2,1),Pos(2,2)}(2);...
	 M{Pos(1,1),Pos(1,2)}{1,1}(1,:),D{Pos(1,1),Pos(1,2)}(1);...
	 M{Pos(1,1),Pos(1,2)}{2,1}(1,:),D{Pos(1,1),Pos(1,2)}(1)];
 B = [M{Pos(3,1),Pos(3,2)}{2,1}(1,:),D{Pos(3,1),Pos(3,2)}(1);...
     M{Pos(3,1),Pos(3,2)}{1,1}(1,:),D{Pos(3,1),Pos(3,2)}(1);...
	 M{Pos(3,1),Pos(3,2)}{2,1}(2,:),D{Pos(3,1),Pos(3,2)}(2);...
	 M{Pos(3,1),Pos(3,2)}{1,1}(2,:),D{Pos(3,1),Pos(3,2)}(2);...
	 M{Pos(1,1),Pos(1,2)}{1,1}(2,:),D{Pos(1,1),Pos(1,2)}(2);...
	 M{Pos(1,1),Pos(1,2)}{2,1}(2,:),D{Pos(1,1),Pos(1,2)}(2)];
 OrigA = A;
 OrigB = B;

% notA enthält alle Punkte, die nicht in A enthalten sind
notA = DMdouble;
for i = 1 : 4
    notA(find(ismember(notA,A(i,:),'rows'),1),:) = [];
end
% notA enthält alle Punkte, die nicht in B enthalten sind
notB = DMdouble;
for i = 1 : 4
    notB(find(ismember(notB,B(i,:),'rows'),1),:) = [];
end
% notAB enthält alle Punkte, die nicht in A oder in B enthalten sind
notAB = DMdouble;
for i = 1 : 6
    notAB(find(ismember(notAB,A(i,:),'rows'),1),:) = [];
    notAB(find(ismember(notAB,B(i,:),'rows'),1),:) = [];
end


%% Optimierung
W3d = zeros(Param_KV(2),Param_KV(2));
betaA = zeros(Param_KV(2),2);
betaB = betaA;
for v = 1 : Param_KV(2)
    if v ~= 1
        B(1:4,1:3) = FNC_Drehung3xyz(B(6,1:3), B(5,1:3), B(1:4,1:3), -betaB(v-1));
    end
    betaB(v) = v/Param_KV(2) * 2* pi;
    B(1:4,1:3) = FNC_Drehung3xyz(B(6,1:3), B(5,1:3), B(1:4,1:3), betaB(v));
    notA = [notAB;B(1:4,:)];
    for w = 1 : Param_KV(2)
        if w ~= 1
            A(1:4,1:3) = FNC_Drehung3xyz(A(6,1:3), A(5,1:3), A(1:4,1:3), -betaA(w-1));
        end
        betaA(w) = w/Param_KV(2) * 2* pi;
        A(1:4,1:3) = FNC_Drehung3xyz(A(6,1:3), A(5,1:3), A(1:4,1:3), betaA(w));
        W3d(w,v) = FNC_Abstandsberechnung(Param_KV(1),A(1:4,1:3),notA(:,1:3),A(1:4,4),notA(:,4),Param_KV(3),Param_KV(4));  % Abstand zum ersten Mutterast
    end
end

if Param_KV(1) == 1
    [x,y]=find(W3d==min(min(W3d)));
elseif Param_KV(1) ==2
    [x,y]=find(W3d==max(max(W3d)));
end
OptBeta = [x,y]*2*pi/Param_KV(2);
OptBetaA = OptBeta(1);  % optimaler Verdrehwinkel für Tochterverzweigung A
OptBetaB = OptBeta(2);  % optimaler Verdrehwinkel für Tochterverzweigung B
AOpt = FNC_Drehung3xyz(OrigA(6,1:3), OrigA(5,1:3), OrigA(1:4,1:3), OptBetaA);  % optimal gedrehte Anschlusspunkte A
BOpt = FNC_Drehung3xyz(OrigB(6,1:3), OrigB(5,1:3), OrigB(1:4,1:3), OptBetaB);  % optimal gedrehte Anschlusspunkte B
DMdouble = [DMdouble(1:size(DMdouble,1)-8,:);[AOpt,OrigA(1:4,4)];[BOpt,OrigB(1:4,4)]];
M_A = {[AOpt(2,:);AOpt(4,:)];[AOpt(1,:);AOpt(3,:)]};
M_B = {[BOpt(2,:);BOpt(4,:)];[BOpt(1,:);BOpt(3,:)]};


end

