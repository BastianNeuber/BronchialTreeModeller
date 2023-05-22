function [P,M_A,M_K,P_ZAM] = FNC_Erzeugung_VE(D,Alpha,tau_M,N,tau_S,mode_S,tau_C,Param_ZAM)
%FNC_Erzeugung_VE
%   Eingabewerte:
%       D        Zeilenvektor mit den Durchmessern der Rohrelemente. Dabei ist der
%                erste Eintrag der Durchmesser des Mutterelementes, die
%                Tochterelemente werden von rechts nach links (pos. x -> neg. x)
%                gezählt
%       Alpha    Zeilenvektor aus den Verzweigungswinkeln. Anordnung der Einträge
%                entsprechend D
%       tau_M    Radienverhältnis von Rohrdurchmesser zu Krümmungsradius
%                Die Parameter stellen Zahlen im Polynom zur Berechnung der
%                Krümmungsradien dar.
%       n        Auflösung: Anzahl der Punkte auf einem Kreis eines jeden
%                Rohrelementes
%       tau_S    Parameter für die Sigmoidfunktion
%       mode_S   Modus der Sigmoidfunktion
%       tau_C    Parameter für die Berechnung der Carinalverrundungsradien
%                R_c
%                  sigma_CR(1) = Parameter für Einfluss der Durchmesser
%                  sigma_CR(2) = Parameter für Einfluss der Winkel
%                  sigma_CR(3) = Parameter für die Breite der Halbellipse

%% Definition von Variablen
n      = N(1);
nRE    = size(D,1);                    % Anzahl der Rohrelemente (inkl. ME)
A      = zeros(n/2,3,2,nRE);           % Abschlusspunkte der Punktreihen => erste Punkte einer Punktreihe, die schon auf dem Mantel der Tochterelemente liegen
B      = zeros(n/2,3,2,nRE-2);         % Endpunkte der Krümmungsbögen aller Punktlinien
C      = zeros(n/2,3,nRE-2);           % Mittelpunkte der Bogenstrecken der Carinalverrundungen
F_te   = zeros(nRE,1);                 % Ergebnisse der Fallunterscheidung für jedes Tochterelement: 0 = norm(B) > norm(U), 1 = norm(B) < norm(U)
F_zf   = zeros(n/2,nRE-2);             % Ergebnisse der Fallunterscheidung für jede Zwischenfläche
Gamma  = zeros(n/2,3,nRE-2);           % Winkel zw. x-Achse und Strecke Mittelpunkt Carinalverrundung - Mittelpunkt Krümmungsbogen
L      = zeros(n/2,3,2,nRE-2);         % Punktreihenlängen
M_C    = zeros(n/2,3,nRE-2);           % Mittelpunkte der Carinalverrundungen
M_K    = zeros(nRE,3);                 % Konstruktionsmittelpunkte (Mittelpunkte der Krümmungsbögen)
M_A    = zeros(nRE,3);                 % Mittelpunkte der Abschlusskreise -> AUSGABEWERT ZUR POSITIONIERUNG DER TOCHTERELEMENTE
nPR    = zeros(nRE-2);                 % Anzahl der Punkte auf einer jeden Punktreihe einer Zwischenfläche
P_pos  = zeros((n/2+1)*(n-1),3);       % Punkte der positiven Außenschale
P_neg = zeros((n/2-1)*(n-1),3);        % Punkte der negativen Außenschale
P_ZAM  = cell(nRE-1,1);                % Variable, die die Punkte des Zentralachsmodells speichert
P_ZF  = zeros(n,3,n/2,2,nRE-2);        % Punkte der Carinalflächen
Phi    = abs(Alpha);                   % pos. Winkel zw. x-Achse und Ende des Krümmungsbogens
Psi    = zeros(n/2,2,nRE-2);           % Auflösung der Punktreihen
%psi    = 0;                           % Restliche Länge zum nächsten Punkt auf dieser Punktreihe in diesem Iterationsschritt
R_C    = zeros(n/2,nRE-2);             % Radien der Carinalverrundung
R_K    = zeros(n/2,2,nRE-2);           % Konstruktionsradien (Radien der Krümmungsbögen)
%theta  = 0;                           % Winkel des Vektors zum nächsten Punkt in einer Punktreihe
U      = zeros(n/2,3,2,nRE-2);         % Berührpunkte zwischen den Punktreihenlinien und der Carinalverrundung


%% Bestimmung von Krümmungsmittelpunkten, Zentralachspunkten, Carinalverrundungsradien und Konstruktionsradien in den Carinalverrundungsflächen
% Schleife über alle Tochterelemente
for w = 2 : nRE
    M_K(w,1) = (tau_M{1}*D(w)^tau_M{2} + tau_M{3}*D(w) + abs(tau_M{4}(w-1)) + tau_M{5}*D(1))*sign(Alpha(w));  % Vorzeichen der Verzweigungswinkel sind umgekehrt gegenüber der Ausarbeitung!!!!
end
% Schleife über alle Zwischenräume
for w = 1 : nRE-2
    for j = 1 : n/2
        R_K(j,1,w) = abs(abs(M_K(w+1,1))+sign(Alpha(w+1))*sin(2*pi*(j-0.5)/n)*D(w+1)/2);
        R_K(j,2,w) = abs(abs(M_K(w+2,1))-sign(Alpha(w+2))*sin(2*pi*(j-0.5)/n)*D(w+2)/2);
        B(j,:,1,w) = [-sign(Alpha(w+1))*cos(Phi(w+1)),0,sin(Phi(w+1))]*R_K(j,1,w) + [0,cos(2*(j-0.5)*pi/n),0]*D(w+1)/2 + M_K(w+1,:);
        B(j,:,2,w) = [-sign(Alpha(w+2))*cos(Phi(w+2)),0,sin(Phi(w+2))]*R_K(j,2,w) + [0,cos(2*(j-0.5)*pi/n),0]*D(w+2)/2 + M_K(w+2,:);
        R_C(j,w) = IFNC_Carinalradius(tau_C,[D(w+1), D(w+2); Alpha(w+1), Alpha(w+2)],j,n);
    end
end
% Übergangseinstellung für die Sigmoidfunktion zwischen Außen- und
% Innenfläche: Der Parameter sigma_S wird so erst in der Mitte (eta > 0)
% wirksam, damit die Zwischen- und Außenflächen zusammentreffen.
if mode_S == 3
    sigma_Sy = 1;
elseif mode_S == 2
    sigma_Sy = 0;
elseif mode_S == 1
    sigma_Sy = 1;
elseif mode_S == 4
    sigma_Sy = 1;
end 

%% Fallunterscheidung: Bestimmung erforderlicher gemoetrischer Kennwerte der Carinalverrundungsbögen:
%   Mittelpunkte (M_c)
%   Übergangspunkte (U)
%   Winkel zw. M_c, M_k und x-Achse (Gamma(j,1/2,i))
%   Öffnungswinkel der Carinalverrundung (Gamma(j,3,i))
%   Längen aller Krümmungsbögen innerhalb einer Punktreihe
for w = 1 : nRE-2
    for j = 1 : n/2
        % Gamma(j,1/2,i) unter der Vorraussetzung, dass die Krümmungsradien nicht begrenzt sind
        Gamma(j,1,w) = acos(((R_K(j,1,w)+R_C(j,w))^2+(abs(M_K(w+1,1))+abs(M_K(w+2,1)))^2-(R_K(j,2,w)+R_C(j,w))^2)/(2*(R_K(j,1,w)+R_C(j,w))*(abs(M_K(w+1,1))+abs(M_K(w+2,1)))));
        Gamma(j,2,w) = acos(((R_K(j,2,w)+R_C(j,w))^2+(abs(M_K(w+1,1))+abs(M_K(w+2,1)))^2-(R_K(j,1,w)+R_C(j,w))^2)/(2*(R_K(j,2,w)+R_C(j,w))*(abs(M_K(w+1,1))+abs(M_K(w+2,1)))));
        % Fall 1: Beide Seiten haben einen Abschnitt II
        if Gamma(j,1,w) <= Phi(w+1) && Gamma(j,2,w) <= Phi(w+2) % Beachte: Phi = Betrag von Alpha!
            F_zf(j,w) = 1;
            M_C(j,:,w) = [-sign(Alpha(w+1))*(R_K(j,1,w)+R_C(j,w)*sign(Alpha(w+1)))*cos(Gamma(j,1,w)),0,(R_K(j,1,w)+R_C(j,w)*sign(Alpha(w+1)))*sin(Gamma(j,1,w))] + M_K(w+1,:);
            U(j,:,1,w) = R_C(j,w)*[cos(Gamma(j,1,w)),0,-sin(Gamma(j,1,w))] + M_C(j,:,w);
            U(j,:,2,w) = R_C(j,w)*[-cos(Gamma(j,2,w)),0,-sin(Gamma(j,2,w))] + M_C(j,:,w);
            L(j,2,1,w) = (Phi(w+1)-Gamma(j,1,w))*R_K(j,1,w);
            L(j,2,2,w) = (Phi(w+2)-Gamma(j,2,w))*R_K(j,2,w);
        % Fälle 2 & 4: Die negative Seite hat keinen Abschnitt II
        elseif Gamma(j,2,w) > Phi(w+2)
            F_zf(j,w) = 2;
            [M_C(j,1,w),M_C(j,3,w)] = FNC_CVCase2367(Alpha(w+2),B(j,:,2,w),R_C(j,w),R_K(j,1,w),M_K(w+1,1),F_zf(j,w));
            Gamma(j,1,w) = atan(M_C(j,3,w)/abs((M_K(w+1,1)-M_C(j,1,w))));
            Gamma(j,2,w) = atan(M_C(j,3,w)/abs((M_K(w+2,1)-M_C(j,1,w))));
            % Fall 4: Beide Seiten haben keinen Abschnitt II
            if Gamma(j,1,w) > Phi(w+1)
                F_zf(j,w) = 4;
                M_C(j,:,w) = FNC_CVCase48([Alpha(w+1);Alpha(w+2)],[B(j,:,1,w);B(j,:,2,w)],R_C(j,w));
                U(j,:,1,w) = R_C(j,w)*[cos(Alpha(w+1)),0,-sin(Alpha(w+1))] + M_C(j,:,w);
                U(j,:,2,w) = R_C(j,w)*[-cos(Alpha(w+2)),0,sin(Alpha(w+2))] + M_C(j,:,w);
                Gamma(j,1,w) = atan(M_C(j,3,w)/(abs(M_K(w+1,1)-M_C(j,1,w))));
                Gamma(j,2,w) = atan(M_C(j,3,w)/(abs(M_K(w+2,1)-M_C(j,1,w))));
            % Fall 2: Nur die positive Seite hat einen Abschnitt II
            else
                U(j,:,1,w) = M_C(j,:,w) + [cos(Gamma(j,1,w)),0,-sin(Gamma(j,1,w))]*R_C(j,w);
                U(j,:,2,w) = M_C(j,:,w) + [-cos(Alpha(w+2)),0,sin(Alpha(w+2))]*R_C(j,w);
                L(j,2,1,w) = (Phi(w+1)-Gamma(j,1,w))*R_K(j,1,w);
            end
        % Fälle 3 & 4: Die positive Seite hat keinen Abschnitt II
        else
            F_zf(j,w) = 3;
            [M_C(j,1,w),M_C(j,3,w)] = FNC_CVCase2367(Alpha(w+1),B(j,:,1,w),R_C(j,w),R_K(j,2,w),M_K(w+2,1),F_zf(j,w));
            Gamma(j,1,w) = atan(M_C(j,3,w)/(abs(M_K(w+1,1)-M_C(j,1,w))));
            Gamma(j,2,w) = atan(M_C(j,3,w)/(abs(M_K(w+2,1)-M_C(j,1,w))));
            % Fall 4: Beide Seiten haben keinen Abschnitt II
            if Gamma(j,2,w) > Phi(w+2)
                F_zf(j,w) = 4;
                M_C(j,:,w) = FNC_CVCase48([Alpha(w+1);Alpha(w+2)],[B(j,:,1,w);B(j,:,2,w)],R_C(j,w));
                U(j,:,1,w) = R_C(j,w)*[cos(Alpha(w+1)),0,-sin(Alpha(w+1))] + M_C(j,:,w);
                U(j,:,2,w) = R_C(j,w)*[-cos(Alpha(w+2)),0,sin(Alpha(w+2))] + M_C(j,:,w);
                Gamma(j,1,w) = atan(M_C(j,3,w)/(abs(M_K(w+1,1)-M_C(j,1,w))));
                Gamma(j,2,w) = atan(M_C(j,3,w)/(abs(M_K(w+2,1)-M_C(j,1,w))));
            % Fall 3: Nur die negative Seite hat einen Abschnitt II
            else
                U(j,:,1,w) = M_C(j,:,w) + [cos(Alpha(w+1)),0,-sin(Alpha(w+1))]*R_C(j,w);
                U(j,:,2,w) = M_C(j,:,w) + [-cos(Gamma(j,2,w)),0,-sin(Gamma(j,2,w))]*R_C(j,w);
                L(j,2,2,w) = (Phi(w+2)-Gamma(j,2,w))*R_K(j,2,w);
            end
        end
        if norm(U(j,:,1,w)-U(j,:,2,w)) < 2*R_C(j,w)
            Gamma(j,3,w) = acos((2*R_C(j,w)^2-norm(U(j,:,1,w)-U(j,:,2,w))^2)/(2*R_C(j,w)^2));
        else
            Gamma(j,3,w) = pi;
        end
        C(j,:,w) = M_C(j,:,w) + R_C(j,w)*[cos(min(Alpha(w+1),Gamma(j,1,w))+Gamma(j,3,w)/2),0,-sin(min(Alpha(w+1),Gamma(j,1,w))+Gamma(j,3,w)/2)];
        L(j,1,1,w) = 1/2*Gamma(j,3,w)*R_C(j,w);
        L(j,1,2,w) = L(j,1,1,w);
    end
    % y-Koordinate der Carinalkurve
    for j = 1 : n/2
        C(j,2,w) = abs(((((n+2)/4)-(j-0.5))/((n+2)/4))^tau_C(3))*sign((j-0.5)-(n+2)/4)*sqrt(1-(norm(C(j,:,w))^2)/(norm(C((n+2)/4,:,w))^2))*D(1)/2;
    end
end

%% Bestimmung der Mittelpunkte aller Tochterelemente
% Schleife über alle Tochterelemente
for w = 2 : nRE
    if w == 2
        if norm(U(ceil(n/4),:,1,1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2) >= norm(B(ceil(n/4),:,1,1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2)  
            M_A(w,:) = U(ceil(n/4),:,1,1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2; %     ML(i,:,1) + [cos(Phi(i)),0,-sin(Alpha(i))]*D(i)/2;  
            F_te(w) = 1;
        else
            M_A(w,:) = B(ceil(n/4),:,1,1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2;
        end
    elseif w == nRE
        if norm(U(ceil(n/4),:,2,nRE-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2) >= norm(B(ceil(n/4),:,2,nRE-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2)
            M_A(w,:) = U(ceil(n/4),:,2,nRE-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2;
            F_te(w) = 1;
        else
            M_A(w,:) = B(ceil(n/4),:,2,nRE-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2;
        end
    else
        if norm(U(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2) >= norm(U(ceil(n/4),:,1,w-1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2) && norm(U(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2) >= norm(B(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2)
            M_A(w,:) = U(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2;
            F_te(w) = 1;
        elseif norm(U(ceil(n/4),:,1,w-1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2) >= norm(U(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2) && norm(U(ceil(n/4),:,1,w-1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2) >= norm(B(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2)
            M_A(w,:) = U(ceil(n/4),:,1,w-1) + [cos(Phi(w)),0,-sin(Alpha(w))]*D(w)/2;
            F_te(w) = 1;
        else
            M_A(w,:) = B(ceil(n/4),:,2,w-2) + [-cos(Phi(w)),0,sin(Alpha(w))]*D(w)/2;
        end
    end
end

%% Bestimmung der Abschlusspunkte der Punktreihen
% Abschlusspunkte A
for w = 2 : nRE
    for j = 1 : n/2
        A(j,:,1,w) = M_A(w,:) + [cos(Alpha(w))*sin(2*(j-0.5)*pi/n),cos(2*(j-0.5)*pi/n),-sin(Alpha(w))*sin(2*(j-0.5)*pi/n)]*D(w)/2;
        A(j,:,2,w) = M_A(w,:) + [-cos(Alpha(w))*sin(2*(j-0.5)*pi/n),cos(2*(j-0.5)*pi/n),sin(Alpha(w))*sin(2*(j-0.5)*pi/n)]*D(w)/2;
    end
end

%% Längen Abschnitt III L(*,3,*,*)
for w = 1 : nRE-2
    for m = 1 : 2
        if F_te(w+m) ~= 0
            for j = 1 : n/2
                % Bei der Berechnung der Länge für den Abschnitt III wird
                % nur die Projektion auf die Bifurkaitonsebene
                % berücksichtigt. U enthält keine Z-Koordinaten, A und B
                % enthalten diese. Daher müssen sie für A und B bei der
                % Berechnung eliminiert werden. Es wird diejenige Strecke
                % ausgewählt, die kleiner ist: L_3 = min{A-B,A-U}.
                L(j,3,m,w) = min(norm([A(j,1,3-m,w+m),0,A(j,3,3-m,w+m)]-U(j,:,m,w)),norm([A(j,1,3-m,w+m),0,A(j,3,3-m,w+m)]-[B(j,1,m,w),0,B(j,3,m,w)]));
            end
        end
    end
end

%% Auflösung Psi und Bestimmung der Oberflächenpunkte aller Zwischenflächen
for w = 1 : nRE-2   % für alle Zwischenflächen
    nPR(w) = ceil(n*sum(sum(L(1,:,:,w)))./(sum(L(1,3,:,w))+Phi(w+1)*R_K(1,1,w)+Phi(w+2)*R_K(1,2,w)));   % in Ausarbeitung p: Anzahl der Punkte für jede halbe Punktreihe
    for m = 1 : 2                                                                                       % für eine (pos./neg.) Hälfte einer Zwischenfläche
        for j = 1 : n/2                                                                                 % für alle Punktreihen 
            Psi(j,m,w) = sum(L(j,:,m,w))/(nPR(w)+1);                                                    % Auflösung in dieser Punktreihenhälfte
            k = 0;                                                                                      % Nummer des Punktes
            L_j = L(j,:,m,w);                                                                           % Länge der aktuellen Punktreihenhälfte und Abstand zur Grenzlinie xi_grenz
            psi = Psi(j,m,w);                                                                           % in Ausarbeitung p_j: Länge bis zum nächsten Punkt ("Längenbudget")
            theta = Phi(w+m);                                                                           % Winkel des Vektors zum nächsten Punkt
            % Schleife übr alle Punkte im Abschnitt 3 (Zylindermantel)
            while L_j(3) > 0
                if L_j(3) > psi
                    k = k + 1;
                    % Oberflächenpunkt in der xz-Ebene zeros(n/2,3,2,nRE)
                    P_ZF(k,:,j,m,w) = [A(j,1,3-m,w+m),0,A(j,3,3-m,w+m)] + [-sin(Alpha(w+m)),0,-cos(Alpha(w+m))]*((k-1)*Psi(j,m,w)+psi);
                    P_ZF(k,2,j,m,w) = sign(C(j,2))*IFNC_Sigmoid(abs(C(j,1)-P_ZF(k,1,j,m,w)),sigma_Sy,[abs(2*C(j,2)),abs(2*A(j,2,3-m,w+m))],abs(A(j,1,3-m,w+m)-C(j,1)),mode_S);
                    L_j(3) = L_j(3) - psi;
                    psi = Psi(j,m,w);
                else
                    psi = psi - L_j(3);
                    L_j(3) = 0;
                end
            end
            % Schleife über alle Punkte im Abschnitt 2 (Krümmungsbogen)
            while L_j(2) > 0
                if L_j(2) > psi
                    k = k + 1;
                    theta = theta - psi/R_K(j,m,w);
                    % Oberflächenpunkt in der xz-Ebene
                    P_ZF(k,:,j,m,w) = R_K(j,m,w)*[-sign(Alpha(w+m))*cos(theta),0,sin(theta)] + M_K(w+m,:);
                    P_ZF(k,2,j,m,w) = sign(C(j,2))*IFNC_Sigmoid(abs(C(j,1)-P_ZF(k,1,j,m,w)),sigma_Sy,[abs(2*C(j,2)),abs(2*A(j,2,3-m,w+m))],abs(A(j,1,3-m,w+m)-C(j,1)),mode_S);
                    L_j(2) = L_j(2) - psi;
                    psi = Psi(j,m,w);
                else
                    psi = psi - L_j(2);
                    L_j(2) = 0;
                end
            end
            % Schleife über alle Punkte im Abschnitt 1 (Carinalverrundung)
            if U(j,1,m,w) - M_C(j,1,w) ~= 0    % Regelfall
                theta = atan((U(j,3,m,w)-M_C(j,3,w))/((U(j,1,m,w)-M_C(j,1,w))));
            else    % Dieser Fall soll nur einen Div/0-Error vermeiden, der Fall sollte eigentlich gar nciht vorkommen können.
                theta = pi/2;
            end
            while L_j(1) > 0
                if L_j(1) > psi
                    k = k + 1;
                    theta = theta + (2*m-3)*psi/R_C(j,w);
                    % Oberflächenpunkte in der xz-Ebene
                    P_ZF(k,:,j,m,w) = (-2*m+3)*R_C(j,w)*[cos(theta),0,sin(theta)] + M_C(j,:,w);
                    P_ZF(k,2,j,m,w) = sign(C(j,2))*IFNC_Sigmoid(abs(C(j,1)-P_ZF(k,1,j,m,w)),sigma_Sy,[abs(2*C(j,2)),abs(2*A(j,2,3-m,w+m))],abs(A(j,1,3-m,w+m)-C(j,1)),mode_S);
                    L_j(1) = L_j(1) - psi;
                    psi = Psi(j,m,w);
                else
                    psi = psi - L_j(1);
                    L_j(1) = 0;
                end
            end
        end
    end
end

%% Erhöhte Auflösung der Carinalkurve
if N(2) > 1
    P_C = zeros((N(2)-1)*n,3);
    P_C(1,:) = C(w,:) - (C(2,:)-C(1,:))*(1/N(2));
    for w = 1 : n/2-1
        for j = 1 : N(2)-1
            P_C((1+w)+(w-1)*(N(2)-1),:) = C(w,:) + (C(w+1,:)-C(w,:))*(j/N(2));
        end
    end
else
    P_C = [];
end

%% Bestimmung der Oberflächenpunkte der positiven Außenfläche
L_Apos = [M_K(2,1)*Alpha(2),(M_A(2,3) - M_K(2,1)*sin(Alpha(2)))/cos(Alpha(2))];    % Länge der Mittellinie
l_Apos = sum(L_Apos)/n;    % Auflösung der rechten Außenfläche
L_v = round((n-1) * L_Apos(1)/(sum(L_Apos)));    % Anzahl der Punkte, die auf Abschnitt II der Außenfläche entfallen
for j = 1 : n/2    % Schleife über alle Punktreihen der rechten Außenfläche
    for k = 1 : L_v    % Schleife über Abschnitt II
        phi_s = k/round(n * L_Apos(1)/(sum(L_Apos))) * Phi(2);
        r_k11 = IFNC_Sigmoid(k,tau_S,[D(1),D(2)],n-1,mode_S);
        P_pos((j-1)*(n+1) + k,1) = M_K(2,1) - sign(Alpha(2))*cos(phi_s)*(abs(M_K(2,1))-sign(Alpha(2))*sin(2*(j-0.5)*pi/n)*r_k11);
        r_k12 = IFNC_Sigmoid(abs(abs(P_pos((j-1)*(n+1) + k,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
            sigma_Sy, ...
            [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(2)], ...
            abs(A(j,1,1,2)-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
            mode_S);
        P_pos((j-1)*(n+1) + k,2) = sin(2*(j-0.5)*pi/n)^2*r_k11*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k12;
        P_pos((j-1)*(n+1) + k,3) = sin(phi_s)*(abs(M_K(2,1))+sign(Alpha(2))*sin(2*(j-0.5)*pi/n)*(-r_k11));
    end
    if L_v == 0
        k = L_v;
        r_k11 = D(1)/2;
    end
    for l = 1 : round((n-1) * L_Apos(2)/(sum(L_Apos)))    % Schleife über Abschnitt III
        % Falls alle Punkte auf Abschnitt III entfallen
        if L_v == 0
            r_k21 = IFNC_Sigmoid(l,tau_S,[D(1),D(2)],n-1,mode_S);
            P_pos((j-1)*(n+1) + l,1) = sin(2*(j-0.5)*pi/n)*D(1)/2 + l*l_Apos*sin(Alpha(2)) - cos(Alpha(2))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
            r_k22 = IFNC_Sigmoid(abs(abs(P_pos((j-1)*(n+1) + l,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
                sigma_Sy, ...
                [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(2)], ...
                abs(A(j,1,1,2)-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
                mode_S);
            P_pos((j-1)*(n+1) + l,2) = sin(2*(j-0.5)*pi/n)^2*r_k21*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k22;
            P_pos((j-1)*(n+1) + l,3) = l*l_Apos*cos(Alpha(2)) - sin(Alpha(2))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
        % Falls nicht alle Punkte auf Abschnitt III entfallen
        else
            r_k21 = IFNC_Sigmoid(k+l,tau_S,[D(1),D(2)],n-1,mode_S);
            P_pos((j-1)*(n+1) + k + l,1) = P_pos((j-1)*(n+1) + k,1) + l*l_Apos*sin(Alpha(2)) - cos(Alpha(2))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
            r_k22 = IFNC_Sigmoid(abs(abs(P_pos((j-1)*(n+1) + k + l,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
                sigma_Sy, ...
                [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(2)], ...
                abs(A(j,1,1,2)-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
                mode_S);
            P_pos((j-1)*(n+1) + k + l,2) = sin(2*(j-0.5)*pi/n)^2*r_k21*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k22;
            P_pos((j-1)*(n+1) + k + l,3) = P_pos((j-1)*(n+1) + k,3) + l*l_Apos*cos(Alpha(2)) + sin(Alpha(2))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
        end        
    end
end 

%% Bestimmung der Oberflächenpunkte der negativen Außenfläche
L_Aneg = [M_K(nRE,1)*Alpha(nRE),(M_A(nRE,3) - M_K(nRE,1)*sin(Alpha(nRE)))/cos(Alpha(nRE))]; %L(n/2,3,2,nRE-2)
l_Aneg = sum(L_Aneg)/n;
L_v = round((n-1) * L_Aneg(1)/(sum(L_Aneg)));
for j = 1 : n/2
    for k = 1 : L_v
        phi_s = k/round(n * L_Aneg(1)/(sum(L_Aneg))) * Phi(nRE);
        r_k11 = IFNC_Sigmoid(k,tau_S,[D(1),D(nRE)],n-1,mode_S);
        P_neg((j-1)*(n+1) + k,1) = M_K(nRE,1) - sign(Alpha(nRE))*cos(phi_s)*(abs(M_K(nRE,1))+sign(Alpha(nRE))*sin(2*(j-0.5)*pi/n)*r_k11);
        r_k12 = IFNC_Sigmoid(abs(abs(P_neg((j-1)*(n+1) + k,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
            sigma_Sy, ...
            [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(nRE)], ...
            abs(abs(A(j,1,2,3))-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
            mode_S);
        P_neg((j-1)*(n+1) + k,2) = sin(2*(j-0.5)*pi/n)^2*r_k11*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k12;
        P_neg((j-1)*(n+1) + k,3) = sin(phi_s)*(abs(M_K(nRE,1))+sign(Alpha(nRE))*sin(2*(j-0.5)*pi/n)*r_k11);
    end
    if L_v == 0
        k = L_v;
        r_k11 = D(1)/2;
    end
    for l = 1 : round((n-1) * L_Aneg(2)/(sum(L_Aneg)))
        % Falls alle Punkte auf den geraden Teil entfallen
        if L_v == 0
            r_k21 = IFNC_Sigmoid(l,tau_S,[D(1),D(nRE)],n-1,mode_S);
            P_neg((j-1)*(n+1) + l,1) = -sin(2*(j-0.5)*pi/n)*D(1)/2 + l*l_Aneg*sin(Alpha(nRE)) + cos(Alpha(nRE))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
            r_k22 = IFNC_Sigmoid(abs(abs(P_neg((j-1)*(n+1) + l,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
                sigma_Sy, ...
                [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(nRE)], ...
                abs(abs(A(j,1,2,3))-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
                mode_S);
            P_neg((j-1)*(n+1) + l,2) = sin(2*(j-0.5)*pi/n)^2*r_k21*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k22;
            P_neg((j-1)*(n+1) + l,3) = l*l_Aneg*cos(Alpha(nRE)) - sin(Alpha(nRE))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
        % Falls nicht alle Punkte auf den geraden Teil entfallen
        else
            r_k21 = IFNC_Sigmoid(k+l,tau_S,[D(1),D(nRE)],n-1,mode_S);
            P_neg((j-1)*(n+1) + k + l,1) = P_neg((j-1)*(n+1) + k,1) + l*l_Aneg*sin(Alpha(nRE)) + cos(Alpha(nRE))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);
            r_k22 = IFNC_Sigmoid(abs(abs(P_neg((j-1)*(n+1) + k + l,1)) - abs(sin(2*(j-0.5)*pi/n)*D(1)/2)), ...
                sigma_Sy, ...
                [cos(2*(j-0.5)*pi/n)*D(1),cos(2*(j-0.5)*pi/n)*D(nRE)], ...
                abs(abs(A(j,1,2,3))-sin(2*(j-0.5)*pi/n)*D(1)/2), ...
                mode_S);
            P_neg((j-1)*(n+1) + k + l,2) = sin(2*(j-0.5)*pi/n)^2*r_k21*cos(2*(j-0.5)*pi/n) + cos(2*(j-0.5)*pi/n)^2*r_k22;
            P_neg((j-1)*(n+1) + k + l,3) = P_neg((j-1)*(n+1) + k,3) + l*l_Aneg*cos(Alpha(nRE)) - sin(Alpha(nRE))*(r_k11-r_k21)*sin(2*(j-0.5)*pi/n);       
        end
    end
end 

%% Zentralachsmodell
for w = 2 : nRE
    % Schleife über alle Punkte eines Astes im Zentralachsmodell:
    % Schrittweite D(w)/M_K(w,1)/Param_ZAM(2) ist der Winkelabschnitt für 
    % ein D Bogenlänge. Wird nur bei Kollisionsdetektion benötigt.
    if Param_ZAM(1) == 1
        % Punkteberechnung der Zentralachse (Abschnitt II)
        for z = 0 : D(w)/M_K(w,1)/Param_ZAM(2) : Alpha(w)
            P_ZAM{w-1}(size(P_ZAM{w-1},1)+1,:) = M_K(w,:) + [-cos(z)*M_K(w,1),0,sin(z)*M_K(w,1)];
        end
        % Letzten Punkt des Kreisbogens einfügen:
        P_ZAM{w-1}(size(P_ZAM{w-1},1)+1,:) = M_K(w,:) + [-cos(Alpha(w))*M_K(w,1),0,sin(Alpha(w))*M_K(w,1)];
        % Punkteberechnung der Zentralachse (Abschnitt III)
%         Halt = 1;
        if F_te(w) == 1
            P_ZAM_III = (0:norm(P_ZAM{w-1}(1,:)-P_ZAM{w-1}(2,:)):norm(P_ZAM{w-1}(size(P_ZAM{w-1},1),:)-M_A(w,:)))';
            P_ZAM_III = [zeros(size(P_ZAM_III,1),2),P_ZAM_III];
            P_ZAM_III = FNC_Drehung3(P_ZAM_III,'y',[0,-Alpha(w),0]) + ones(size(P_ZAM_III,1),1)*P_ZAM{w-1}(size(P_ZAM{w-1},1),:);
            P_ZAM{w-1} = [P_ZAM{w-1};P_ZAM_III];
        end
        % Durchmesserberechnung zu jedem Punkt der Zentralachse
    end
end


%% Ausgabe
P = P_pos;
for w = 1 : nRE-2
    for j = 1 : n/2
        for m = 1 : 2
            P = [P;P_ZF(:,:,j,m,w)];
        end
    end
end
P = [P;P_neg;P_C];
% Nullzeilen löschen
P(sum(P~=0,2)==0,:) = [];

%% Testplots
% %_________________________________________________________________________________________________________________________________________________
%     close all
%     hold on
%     axis equal
% %     scatter3(P_neg(:,1),P_neg(:,2),P_neg(:,3),2);
%     scatter3(P_pos(:,1),P_pos(:,2),P_pos(:,3),2);
%     scatter3(P(5253:end-5354,1),P(5253:end-5354,2),P(5253:end-5354,3),2);
%     scatter3(P(:,1),P(:,2),P(:,3),1);    
%     p = plot3(C(:,1),C(:,2),C(:,3));
%     p.Color = [0,0,0];
%     for o = 1 : 2
%         for mm = 1 : 2
%             scatter3(A(:,1,mm,o+1),A(:,2,mm,o+1),A(:,3,mm,o+1),2,'r');
%         end
%     end
%     view([90,0,0]);
%     axis equal
% P = FNC_Drehung3(P,'y',[0,-pi/2,0]);
%     Halt = 1;
% %_________________________________________________________________________________________________________________________________________________
end