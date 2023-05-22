function [P] = FNC_Zentralachsmodell(mode,L,D,Param,Alpha,P_VE,Beta)

% Zentralachse eines Zylinders
if mode == 1
    P = 0 : D/Param : L;
    P = P';
    P = [zeros(size(P,1),2),P,ones(size(P,1),1)*D];
        
% Zentralachsen für mehrere Zylinder an einem Verzweigungselement
elseif mode == 2
    P = [];
    I = size(P_VE,1);
    for i = 1 : I
        p = 0 : D(i+1)/Param : L(i)-0.000001;
        p = p';
        p = [zeros(size(p,1),2),p,ones(size(p,1),1)*D(i+1)];
        P_VE{i+I} = [FNC_Drehung3xyz([0,0,0], [0,1,0], p(:,1:3), Alpha(i)) + ones(size(p,1),1)*P_VE{i}(size(P_VE{i},1),:),p(:,4)];
        P = [P;[P_VE{i},ones(size(P_VE{i},1),1)*D(1)];P_VE{i+I}];
    end
end
P = [FNC_Drehung3xyz([0,0,0],[0,0,-1],P(:,1:3),Beta),P(:,4)];

end

