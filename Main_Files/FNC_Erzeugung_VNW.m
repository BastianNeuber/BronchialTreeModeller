function [P,Errors] = FNC_Erzeugung_VNW(N,Gen,ME,Param_RE,Param_VE,Param_Modell,Eingabemodus,Param_PP,Param_WS,Param_KV,M_soll,Param_AE)
% Die Funktion ruft die verschiedenen weiteren Funktion zur Erzeugung des
% Rohrnetzwerkes auf.

% Argumente:
%   n               Auflösung des Modells (Punkte pro KReisquerschnitt)
%   Gen             Anzahl der Generationen, die berechnet werden sollen. Dabei
%                   zählt das erste Mutterelement noch nicht als Generation mit
%   ME              Geometrieparameter des ersten Mutterlementes
%                       ME(1)   Durchmesser des ersten Mutterelementes
%                       ME(2)   Länge des ersten Mutterelementes
%   Param_RE        Cell-Variable mit geometrischen Eigenschaften der
%                   Rohrelemente:
%                      Param_RE{1} = lim_nTE       Grenzen (min,max) der Anzahl von
%                                                  Tochterelementen eines
%                                                  Mutterelementes
%                       Param_RE{2} = Param_D      Parameter der Berechnung von
%                                                  Durchmessern der Tochterelemente
%                                                    Param_D(1) = Varianz
%                                                    Param_D(2) = Erweiterungsexponent
%                       Param_RE{3} = Param_L      Parameter der Berechnung von
%                                                  Längen der Tochterelemente
%                                                    Param_L(1) = Varianz
%                                                    Param_L(2) = Erwartungswert
%                       Param_RE{4} = Param_Alpha  Grenzen (min,max,nenn) der Verzweigungswinkel
%                       Param_RE{5} = VarBL        Länge einer Bronchie und
%                                                  Längenverteilung
%                                                    Param_RE{5}(1) = Ordnung
%                                                    Param_RE{5}(2) = Faktor
%                                                    Param_RE{5}(3) = Anteil
%   Param_VE        Cell-Variable mit geometrischen Eigenschaften der
%                   Verzweigungselemente:
%                       Param_VE{1} = Param_K      Verhältnis von Krümmungsradius
%                                                  zum Rohrdurchmesser
%                       Param_VE{2} = Param_C      Parameter für die
%                                                  Carinalverrundung
%                       Param_VE{3} = Param_S      Parameter für die Sigmoidfunktion:
%                                                    Param_S(1) = Modus:
%                                                       (1 = Tangens Hyperbolicus, 2 = Methode nach Lee et al.)
%                                                    Param_S(2) = Modulationsparameter
%   Param_Modell    Parameter für das Gesamtmodell
%                       Param_Modell{1}            Dimension des Gesamtmodells
%                                                    Param_Modell{1} = 2: Alle Elemente liegen in der Bifurkationsebene
%                                                    Param_Modell{1} = 3: Der Bronchialbaum ragt aus der Bifurkationsebene hinaus => notwendig zur Kollisionsvermeideung
%   Eingabemodus    Speichert, ob vorgefertigte Eingabeparameter für die
%                   geometischen Eigenschaften der Rohrelemente verwendet werden ('m'), oder ob
%                   diese Werte zufällig bestimmt werden ('a'). Für den
%                   Durchlauf zur Berechnung der Wandstärke wird der
%                   Eingabemodus 'w' verwendet.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Erster Teil: Initialisierung und Modellvorbereitung  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all
clc

%% Initialisierung von Variablen
if Eingabemodus == 'a'
    D = cell(Gen,1);
    L = D;
    Alpha = D;
    Beta = D;
    D{1,1} = ME(1);
    L{1,1} = ME(2);
    Alpha{1,1} = [0,0,0];
    M_K = cell(Gen,1);
    BL = D;
elseif Eingabemodus == 'm'
    D = Param_RE{1};
    L = Param_RE{2};
    % Falls die Bronchienlängen eingehalten werden sollen, stellt die
    % Matrix L die Bronchienlänge statt der Länge der Rohrelemente dar. In
    % diesem Fall muss sie als Bronchienlänge umdefiniert werden, damit der
    % Code allgemeingültig bleibt.
    if Param_VE{1}{6} == true
        BL = L;
    end
    Alpha = Param_RE{3};
    Beta = Param_RE{4};
    Gen = size(D,1);
    M_K = cell(Gen,1);
elseif Eingabemodus == 'w'
    D = Param_RE{1};
    D{1,1} = D{1,1} + 2*Param_WS(2)*Param_WS(3);
    L = Param_RE{2};
    Alpha = Param_RE{3};
    Beta = Param_RE{4};
    Gen = size(D,1);
    M_K_WS = Param_VE{1};
    Param_VE2_nenn = Param_VE{2};
end
Errors = false(999,1);
n_VZW = 1;
M = cell(Gen,1);
M{1,1} = cell(2,1);
M{1,1}{1} = [0,0,0];
M{1,1}{2} = [0,0,ME(2)];
M_RE = cell(2,1);
Koll = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Zweiter Teil: Erzeugung der Hauptgeometrieparameter  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for g = 2 : Gen  % Schleife über alle Generationen
    Offset = 0;
    for v = 1 : n_VZW  % Schleife über alle Verzweigungen einer Generation
        for t = 1 : size(D{g-1,v},1)  % Schleife über alle Tochteräste einer Verzweigung
            %% Anpassung von Geometrieparametern
            if Eingabemodus == 'a'  % Falls der stochastische Erzeugungsprozess zur Parameterwerterzeugung verwendet wird, ...
                [D{g,Offset + t},L{g,Offset + t},Alpha{g,Offset + t},Beta{g,Offset + t}] = FNC_Erzeugung_Parameter(D{g-1,v}(t),Param_RE{1},Param_RE{4},Param_RE{2},Param_RE{3},Param_Modell(1));
                Alpha{g,Offset + t} = [zeros(size(Alpha{g,Offset + t},1),1),Alpha{g,Offset + t},zeros(size(Alpha{g,Offset + t},1),1)];
            elseif  Eingabemodus == 'm'  % Falls der manuelle Eingabemodus verwendet wird und die Brochianlängen eingehalten werden sollen...
                if Param_VE{1}{6} == true
                    Param_VE{1} = {0,0,0,Param_VE{1,1}{7}*BL{g,Offset + t}./Alpha{g,Offset + t}(:,2),0};
                    L{g,Offset + t} = Param_VE{1,1}{7}*BL{g,Offset + t};
                end
            end
            %% Prüfung auf Fehler EC003
            if Param_WS(3) == -1 && D{g,Offset + t}(1) < 2*Param_WS(2) && D{g,Offset + t}(2) < 2*Param_WS(2)
                Errors(3) = true;
                P = [];
                save('VAR_temp.mat');
                APP_Ergebnis_speichern
                return
            end
        end
        %% Ermitteln, an welcher Stelle in der Variable D gelesen werden muss
        Offset = Offset + size(D{g-1,v},1);  % Diese Zeile verallgemeinert den dichotomen Verzweigungsfall und ermittelt die aktuelle Spalte in der Variablen D.
    end
    %% Anzahl der Verzweigungen in der aktuellen Generation bestimmen
    empty = cellfun('isempty',D);
    n_VZW = sum(empty(g,:)==0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Dritter Teil: Kollisionsvermeidung  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Param_KV(1) ~= 0 || Param_KV(5) == 1) && Eingabemodus ~= 'w'  % Falls alle Bronchien in die Kollisionsvermeidung einbezogen werden sollen...
    DMdouble = [0,0,0,D{1,1};0,0,L{1},D{1,1}];  % Durchmesser und Anschlussmittelpunkte als double-Variable

%     %% Erstes Mutterelement erzeugen
%     P = FNC_Erzeugung_RE(N(1),ME(2),D{1,1});

    %% Bronchialbaummodell erzeugen
    for g = 2 : Gen  % Schleife über alle Generationen
        Offset = 0;
        for v = 1 : n_VZW  % Schleife über alle Verzweigungen einer Generation
            for t = 1 : size(D{g-1,v},1)  % Schleife über alle Tochteräste einer Verzweigung

                %% Erzeugung der Anschlusspunkte
                M_RE{1} = FNC_Erzeugung_VE_Anschlusspunkte([D{g-1,v}(t);D{g,Offset + t}],[0;Alpha{g,Offset + t}(:,2)],Param_VE{1},[N(1),N(2)],Param_VE{3}(2),Param_VE{3}(1),Param_VE{2});
                M_RE{1} = M_RE{1}(2:end,:);  % Anpassung des Variablenformats: Löschen der ersten Zeile, da FNC_Erzeugung_VE eine Nullzeile für den Mutterast produziert.

                %% Bestimmung der lokalen Anschlussmittelpunkte
                for tt = 1 : size(M_RE{1},1)  % M_RE{1} ist schon bekannt, M_RE{2} muss aus alpha und Länge des Rohrelementes noch berechnet werden.
                    M_RE{2}(tt,:) = M_RE{1}(tt,:) + [sin(Alpha{g,Offset + t}(tt,2)),0,cos(Alpha{g,Offset + t}(tt,2))]*L{g,Offset + t}(tt);
                end

                %% Platzieren der Anschlusspunkte im globalen Koordinatensystem
                M{g,Offset + t} = cell(2,1);
                for s = 1 : 2
                    M{g,Offset + t}{s} = FNC_PlatzierenGruppe(M{g-1,v}{1}(t,:),M{g-1,v}{2}(t,:),M_RE{s});
                end
                DMdouble = [DMdouble;M{g,Offset + t}{1},D{g,Offset + t};M{g,Offset + t}{2},D{g,Offset + t}];

                %% Optimierung des Verdrehwinkels
                % Falls...
                %   ...mindestens die dritte Generation berechnet wird,
                %   ...eine gerade Anzahl von Verzweigungen in dieser Generation berechnet wurde
                if g >= 3 && rem(t,2) == 0 && Param_KV(1) ~= 0
                    [Beta{g,Offset + 1},Beta{g,Offset + 2},DMdouble,M{g,Offset + 1},M{g,Offset + 2}] = FNC_OptBeta(D,M,DMdouble,Param_KV);
                end
                
            end
            %% Ermitteln, an welcher Stelle in der Variable D gelesen werden muss
            Offset = Offset + size(D{g-1,v},1);
        end
        %% Anzahl der Verzweigungen in der aktuellen Generation bestimmen
        empty = cellfun('isempty',D);
        n_VZW = sum(empty(g,:)==0);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Vierter Teil: Oberflächenpunkte  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DMdouble = [0,0,0,D{1,1};0,0,L{1},D{1,1}];  % Durchmesser und Anschlussmittelpunkte als double-Variable

%% Erstes Mutterelement erzeugen
P = FNC_Erzeugung_RE(N(1),ME(2),D{1,1});

%% Zentralachsmodell des ersten Mutterelementes erzeugen
if Param_KV(5) == 1
    P_ZAM = FNC_Zentralachsmodell(1,ME(2),D{1,1},Param_KV(6),[],[],0);
end

%% Bronchialbaummodell erzeugen
for g = 2 : Gen  % Schleife über alle Generationen
    Offset = 0;
    for v = 1 : n_VZW  % Schleife über alle Verzweigungen einer Generation
        for t = 1 : size(D{g-1,v},1)  % Schleife über alle Tochteräste einer Verzweigung
            
            %% Anpassung der Konstruktionsmittelpunkte
            %  (nur bei der Konstruktion der Wandstärke)
            if Eingabemodus == 'w'  % Falls die Wandstärke konstruiert wird, ...
                D_nenn = D{g,Offset + t};
                D{g,Offset + t} = D{g,Offset + t} + ones(size(D{g,Offset + t},1),1)*2*Param_WS(2)*Param_WS(3);  % aktuelle Durchmesser aktualisieren
                Param_VE{1} = {0,0,0,M_K_WS{g,Offset + t}(2:end,1),0};  % Konstruktionsmittelpunkte fixieren (nur konstanter Anteil)
                % Carinalverrundungsparameter sowohl für den Durchmesser als auch für den Verzweigungswinkel umgekehrt proportional von der Durchmesseränderung 
                % des Mutterastes durch Wandstärke abhängig machen. => Außenwand führt zu engerer Carinalverrundung, Innenwand zu größerer Carinalverrundung
                if Param_WS(4) == 1
                    Param_VE{2} = [sum(D_nenn)^2/(sum(D{g,Offset + t})^2)*Param_VE2_nenn(1), ...  % Abhängigkeit vom Durchmesser
                        sum(D_nenn)^2/(sum(D{g,Offset + t})^2)*Param_VE2_nenn(2), ...  % Abhängigkeit vom Verzweigungswinkel
                        Param_VE2_nenn(3)];  % Breite der Carinalkurve bleibt unverändert.
                end
            end
            
            %% Erzeugung des Verzweigungselementes
            [P_VZW,M_RE{1},M_K{g,Offset + t},p_ZAM] = FNC_Erzeugung_VE([D{g-1,v}(t);D{g,Offset + t}],[0;Alpha{g,Offset + t}(:,2)],Param_VE{1},[N(1),N(2)],Param_VE{3}(2),Param_VE{3}(1),Param_VE{2},[Param_KV(5),Param_KV(6)]);

            %% Vervollständigung des lokalen Zentralachsmodells
            %  Das ZAM wurde im vorhergegangenen Schritt nur für das VE und
            %  nicht für die RE konstruiert. Das zusammengesetzte ZAM wird
            %  anschließend mit dem Verdrehwinkel gedreht.
            if Param_KV(5) == 1 && Eingabemodus ~= 'w'
                p_ZAM = FNC_Zentralachsmodell(2,L{g,Offset + t},[D{g-1,v}(t);D{g,Offset + t}],Param_KV(6),Alpha{g,Offset + t}(:,2),p_ZAM,Beta{g,Offset + t});
            end
            M_RE{1} = M_RE{1}(2:end,:);  % Anpassung des Variablenformats: Löschen der ersten Zeile, da FNC_Erzeugung_VE eine Nullzeile für den Mutterast produziert.
                
            %% Erzeugung der Rohrelemente
            for r = 1 : size(D{g,Offset + t},1)
                P_RE{r} = FNC_Erzeugung_RE(N(1),L{g,Offset + t}(r),D{g,Offset + t}(r));
            end
                            
            %% Bestimmung der lokalen Anschlussmittelpunkte
            for tt = 1 : size(M_RE{1},1)  % M_RE{1} ist schon bekannt, M_RE{2} muss aus alpha und Länge des Rohrelementes ncoh berechnet werden.
                M_RE{2}(tt,:) = M_RE{1}(tt,:) + [sin(Alpha{g,Offset + t}(tt,2)),0,cos(Alpha{g,Offset + t}(tt,2))]*L{g,Offset + t}(tt);
            end

            %% Platzierung der Rohrelemente im lokalen Koordinatensystem und Drehung der Verzweigung im lokalen Koordinatensystem
            [P_Gruppe,M_REx] = FNC_PlatzierenRE([M_RE{1};M_RE{2}],Alpha{g,Offset + t},P_VZW,P_RE,Beta{g,Offset + t});
            M_RE{1} = M_REx(1:2,:);  % M_RE aktualisieren (M_REx enthält Verdrehwinkel jeder Verzweigung)
            M_RE{2} = M_REx(3:4,:);  % M_RE aktualisieren (M_REx enthält Verdrehwinkel jeder Verzweigung)

            %% Platzieren der Gruppe (Oberflächenpunkte) im globalen Koordinatensystem
            P = [P;FNC_PlatzierenGruppe(M{g-1,v}{1}(t,:),M{g-1,v}{2}(t,:),P_Gruppe)];
            
            %% Platzieren der Gruppe (Zentralachsmodell) im globalen Koordinatensystem und Drehung mit Verdrehwinkel
            if Param_KV(5) == 1
                size_P_ZAM_alt = size(P_ZAM,1);
                P_ZAM = [P_ZAM;[FNC_PlatzierenGruppe(M{g-1,v}{1}(t,:),M{g-1,v}{2}(t,:),p_ZAM(:,1:3)),p_ZAM(:,4)]]; % [p_ZAM{1};p_ZAM{2};p_ZAM{3};p_ZAM{4}])];
                % Richtungsvektor
                Richtungsvektor = FNC_Richtungsvektor(M{g-1,v}{1}(t,:),M{g-1,v}{2}(t,:));
            end

            %% Platzieren der Anschlusspunkte im globalen Koordinatensystem
            M{g,Offset + t} = cell(2,1);
            for s = 1 : 2
                M{g,Offset + t}{s} = FNC_PlatzierenGruppe(M{g-1,v}{1}(t,:),M{g-1,v}{2}(t,:),M_RE{s});
            end
            DMdouble = [DMdouble;M{g,Offset + t}{1},D{g,Offset + t};M{g,Offset + t}{2},D{g,Offset + t}];
            
            %% Prüfung auf den Fehler EC001
            %  (Anschlusspunktverschiebung bei Wandstärkekonstruktion)
            if Eingabemodus == 'w'
                if norm(M{g,Offset + t}{1} - M_soll{g,Offset + t}{1}) > Param_WS(5) || norm(M{g,Offset + t}{2} - M_soll{g,Offset + t}{2}) > Param_WS(5)
                    Errors(1) = true;
                end
            end
            
            %% Prüfung auf den Fehler EC002
            %  (Kollision)
            if Param_KV(5) == 1 && g >= 3 && Eingabemodus ~= 'w'
                Koll = FNC_Kollisionspruefung(P_ZAM,Param_WS(1:3),Param_KV(7),size_P_ZAM_alt,Koll,Richtungsvektor);
            end
            
        end
        %% Ermitteln, an welcher Stelle in der Variable D gelesen werden muss
        Offset = Offset + size(D{g-1,v},1);
    end
    %% Anzahl der Verzweigungen in der aktuellen Generation bestimmen
    empty = cellfun('isempty',D);
    n_VZW = sum(empty(g,:)==0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Fünfter Teil: Nachbereitung  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Post-Processing
if Eingabemodus ~= 'w'
    
    %% Wandstärke
    if Param_WS(1) == 1
        [P_WS,Errors] = FNC_Erzeugung_VNW(N,...  % Nennauflösung
            Gen,...  % Anzahl Generationen
            ME,...  % Durchmesser und Länge des ersten Mutterastes
            {D,L,Alpha,Beta},...  % Hauptgeometrieparameter
            {M_K,Param_VE{2},Param_VE{3}},...  % Parameter für die Verzweigungselemente
            0,...  % Einstellung zu Bifurkationsebenen
            'w',...  % Eingabemodus
            Param_PP,...  % Parameter zum Post Processing
            [0,Param_WS(2),Param_WS(3),Param_WS(4),Param_WS(5)],...  % Parameter zur Wandstärke
            [0,0,0,0,0,0],...  % Einstellungen zur Kollisionsvermeidung
            M,...  % Soll-Anschlussmittelpunkte zur Erkennung von Fehlern bei der Wandstärkekonstruktion
            0);  % Parameter zu Endstücken
        P = [P;P_WS];
    end
    
    %% Sockel
    if Param_PP{3} == 1  % Falls ein Sockel erzeugt wird...
        P = [P;FNC_Erzeugung_Sockel(max(D{1,1}+2*Param_WS(2)*Param_WS(1)*Param_WS(3),Param_WS(1)*D{1,1}), ...  % D0 (Mutterast außen)
                Param_PP{4}, ...  % L (Sockel)
                Param_PP{5}, ...  % D (Sockel)
                N(3), ...         % n (Sockelauflösung)
                min(D{1,1}+2*Param_WS(2)*Param_WS(1)*Param_WS(3),Param_WS(1)*D{1,1}))];  % B = Bohrung (Mutterast innen)
    elseif Param_PP{3} == 0 && Param_WS(1) ~= 0 && Param_PP{6} == 1  % Falls kein Sockel, aber eine Wandstärke erzeugt wird, muss ein unterer Deckel berechnet werden.
        P = [P;FNC_Erzeugung_Deckel(N(1), ...  % n (Auflösung)
                min(D{1,1},D{1,1} + Param_WS(2)*2*Param_WS(1)*Param_WS(3)), ...  % Di (Innendurchmesser)
                max(D{1,1},D{1,1} + Param_WS(2)*2*Param_WS(1)*Param_WS(3)))];  % Da (Außendurchmesser)
    elseif Param_PP{3} == 0 && Param_WS(1) == 0 && Param_PP{6} == 1
        P = [P;FNC_Erzeugung_Deckel(N(1),0,D{1,1})];
    end

    %% Endstücke und Deckel
    if Param_PP{1} == 1  % Option Button "gerade Zylinder" im Postprocessingfeld ist aktiv
        for i = 1 : 2^(Gen-1)
            P_AE = FNC_Erzeugung_ESzyl(N(1), ...  %Auflösung
                Param_AE{2}, ...  % Länge der Abschlusszylinder
                min(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)), ...  % Di (Innendurchmesser)
                max(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)), ...  % Da (Außendurchmesser)
                Param_PP{6});  % soll das Modell geschlossen werden?
            P = [P;FNC_PlatzierenGruppe(M{g,ceil(i/2)}{1}(2-mod(i,2),:),M{g,ceil(i/2)}{2}(2-mod(i,2),:),P_AE)];
        end
    elseif Param_PP{1} == 0 && Param_PP{6} == 1 && Param_WS(1) == 0  % Option Button "keine" ist aktiv, das Modell soll geschlossen werden und es gibt keine Wandstärke
        for i = 1 : 2^(Gen-1)
            P_Deckel = FNC_Erzeugung_Deckel(N(1),0,D{g,ceil(i/2)}(2-mod(i,2)));
            P = [P;FNC_PlatzierenGruppe(M{g,ceil(i/2)}{1}(2-mod(i,2),:),M{g,ceil(i/2)}{2}(2-mod(i,2),:),P_Deckel)];
        end
    elseif Param_PP{1} == 0 && Param_PP{6} == 1 && Param_WS(1) == 1  % Option Button "keine" ist aktiv, das Modell soll geschlossen werden und es gibt eine Wandstärke
        for i = 1 : 2^(Gen-1)
            P_Deckel = FNC_Erzeugung_Deckel(N(1), ...
                min(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)), ...  % Di (Innendurchmesser)
                max(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)));     % Da (Außendurchmesser)
            P = [P;FNC_PlatzierenGruppe(M{g,ceil(i/2)}{1}(2-mod(i,2),:),M{g,ceil(i/2)}{2}(2-mod(i,2),:),P_Deckel)];
        end
    elseif Param_PP{1} == 2  % Option Button "benutzerdefiniert" ist aktiv
        [P_ESi,P_ESa,P_ESd] = FNC_Erzeugung_ESbd3D(Param_AE{2},Param_AE{3},Param_AE{4});  % Funktion erstellt dreidimensionales Endstück (aus drei Sets von Punkten: Innen, Außen, Deckel)). Dabei wird zunächst der veränderliche Durchmesser der Terminalbronchien nicht berücksichtigt (=> vollständiger Kegel).
        for i = 1 : 2^(Gen-1)
            % Funktion zur Anpassung auf variable Durchmesser und Formatierung der Oberflächenpunkte der Endstücke
            P_Endstueck = FNC_Anpassung_bdES(N(1),...  % Auflösung der Bohrung
                min(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)),...  % Di (Innendurchmesser)
                max(D{g,ceil(i/2)}(2-mod(i,2)),D{g,ceil(i/2)}(2-mod(i,2)) + Param_WS(2)*2*Param_WS(1)*Param_WS(3)),...  % Da (Außendurchmesser)
                P_ESi,P_ESa,P_ESd);  % Punkte des allgemeingültigen dreidimensionalen Endstückes
            if isempty(P_Endstueck) == false
                P = [P;FNC_PlatzierenGruppe(M{g,ceil(i/2)}{1}(2-mod(i,2),:),M{g,ceil(i/2)}{2}(2-mod(i,2),:),P_Endstueck)];
            end
        end
    end
    

%% Prüfung auf Fehler EC002
%     hold off
    if isempty(Koll) == false
        Errors(2) = true;
        hold on
        scatter3(Koll(:,1),Koll(:,2),Koll(:,3),200,[1,0,0],'filled');
    end
    
%% Visualisierung und Workspace speichern
%     figure
    scatter3(P(:,1),P(:,2),P(:,3),0.7,[0, 0.4470, 0.7410]);
    axis equal
    view([1,1,0.5]);
% % Einblenden ZAM (Anfang)
%     hold on
%     scatter3(P_ZAM(:,1),P_ZAM(:,2),P_ZAM(:,3),100,[0.6,0.2,0.4]);
% % Einblenden ZAM (Ende)

%% Ergebnisse speichern
    save('VAR_temp.mat');
    APP_Ergebnis_speichern
end
