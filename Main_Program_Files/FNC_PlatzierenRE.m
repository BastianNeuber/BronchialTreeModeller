function [P,M] = FNC_PlatzierenRE(M,Alpha,P_VZW,P_RE,beta)
% Die Funktion verschiebt alle Rohrelemente zu einem Verzweigungselement an
% die richtige Position und schreibt alle Punkte in eine gemeinsame Matrix.

P = P_VZW;
if isempty(P_RE) == 0  % Falls nur Anschlusspunkte berechnet werden sollen, wird von der Hauptfunktion eine leere Matrix übergeben.
                       % Die Transformation der Oberflächenpuntke wird dann übersprungen.
    for i = 1 : size(P_RE,2)
        % Drehung der RE um die lokale z-Achse mit dem Verzweigungswinkel
        Q = FNC_Drehung3(P_RE{i},'y',-Alpha(i,:));
        % Verschiebung der RE an die Mittelpunkte der Abschlusskreise des Verzweigungselementes
        Q = Q + ones(size(Q,1),1)*M(i,:);
        P = [P;Q];
    end
    P = FNC_Drehung3(P(:,:),'z',[0,0,beta]);
end

M = FNC_Drehung3(M,'z',[0,0,beta]);

end