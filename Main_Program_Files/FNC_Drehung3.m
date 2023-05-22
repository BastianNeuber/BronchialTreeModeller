function [P] = FNC_Drehung3(P,O,phi)
%DREHUNG_Z dreht die eingegebene Punktewolke P um die Achse O um den
%eingegebenen Winkel phi

%% Sicherung für den zweidimensionalen Fall

if O == 'x' || O == 'a'
    R = [1, 0, 0; 0, cos(phi(1)), -sin(phi(1)); 0, sin(phi(1)), cos(phi(1))];
    P = P*R;
end
if O == 'y' || O == 'a'
    R = [cos(phi(2)), 0, sin(phi(2)); 0, 1, 0; -sin(phi(2)), 0, cos(phi(2))];
    P = P*R;
end
if O == 'z' || O == 'a'
    R = [cos(phi(3)), -sin(phi(3)), 0; sin(phi(3)), cos(phi(3)), 0; 0, 0, 1];
    P = P*R;
end

end