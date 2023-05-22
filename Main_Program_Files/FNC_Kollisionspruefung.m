function [Koll] = FNC_Kollisionspruefung(P,WS,Puffer,size_alt,Koll,RV)

% P        Punkte des ZAM
% WS       Wandstärke
%   WS(1)  Existenz der WS
%   WS(2)  Betrag der WS
%   WS(3)  Richtung der WS
% Puffer   Sicherheitsabstand zwischen Bronchien
% size_alt Zeilennummer, bis zu der die Punkte in P zum alten, also nicht
%          aktuellen ZAM gehören (die aktuelle Verzweigung istr diejenige,
%          mit der die Kollisionen detektiert werden)
% Koll     Matrix aller kollidierenden Punkte im ZAM
% RV       Richtungsvektor der aktuellen Verzweigung (zur Identifikation
%          irrelevanter Punkte des ZAM, bei denen eine Kollision
%          ausgeschlossen werden kann


% Halbierung des Richtungsparameters der Wandstärke, um Einfluss bei WS
% nach außen zu 1 zu setzen und Einfluss der WS nach innen zu 0 zu setzen.
WS(3) = WS(3)/2;

% scatter3(P(1:size_alt,1),P(1:size_alt,2),P(1:size_alt,3));
% hold on
% scatter3(P(size_alt+1:end,1),P(size_alt+1:end,2),P(size_alt+1:end,3));
% % view([0,-1,0]);
% axis equal

% Schleife über die Punkte des fixen ZAM
for i = 1 : size_alt
    % Winkel zwischen dem Richtungsvektor RV und dem Ortsvektor des zu
    % prüfenden Punktes: Ist der Winkel zwischen diesen beiden größer
    % als 90°, dann kann eine Kollision ausgeschlossen werden: der
    % Punkt liegt dann "hinter" der aktuellen Verzweigung. Der
    % Richtungsvektor des zu untersuchenden Punktes wird vom
    % Verzweigungspunkt (Anschlusspunkt) zur Mutterast der aktuellen
    % Verzweigung berechnet:
    x = P(i,1:3) - P(size_alt+1,1:3);
    phi = acos(dot(x,RV)/(norm(x)*norm(RV)));
    if phi < pi/2
        % Schleife über die Punkte der aktuellen Verzweigung des ZAM
        for j = size_alt+1 : size(P,1)
            %  Abstand der ZAM-Punkte  -    Radiensumme    -              WS             - Puffer
            if norm(P(j,1:3)-P(i,1:3)) - (P(j,4)+P(i,4))/2 - 2*(WS(1)*WS(2)*(WS(3)+0.5)) - Puffer < 0
                Koll = [Koll;P(j,1:3);P(i,1:3)];
            end
        end
    end
end
















