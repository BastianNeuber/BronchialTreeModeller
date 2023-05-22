function [sigma] = IFNC_Sigmoid(x,tau,D,x_max,mode)
%Sigmoidmfunktion zur Bestimmung der Annäherung zweier Durchmesser in der
%Bifurkation. Die Funktion berechnet den Betrag, der zur Berechnung des
%Konstruktionsradius einer Punktreihe vom mittleren Radius (also vom
%x-Achsenabschnitt des Konstruktionsmittelpunktes M_k) subtrahiert bzw.
%addiert werden muss.
%
%sigma(0) = D(1)/2
%sigma(Phi) = D(2)/2
%
%Zwischen diesen beiden Grenzwerten nähert sich der Funktionswert den
%Grenzwerten stetig an. Die Funktion hat genau in der Mitte zwsichen den
%Grenzwerten einen Wendepunkt.
%
%Eingangsparameter:
%   x      Funktionsvariable: Laufvariable der Funktion, die diese Funktion aufruft [rad]
%   beta   Modulator für die Geschwindigkeit des Übergangs [-]
%   D      Vektor aus den Durchmessern der zu verbindenden Rohrelemente:
%            D(1) = Durchmesser des Startelementes
%            D(2) = Durchmesser des Zielelementes
%   x_max  Obere Grenze des Definitionsbereichs der Funktion. Die untere
%          Grenze liegt automatisch immer bei 0.
%   mode   Modus der Sigmoidalfunktion:

if mode == 1  % Im Hauptmenü nicht wählbar/enthalten.
    % Sigmoidfunktion als Tangens-Hyperbolicus
    deltaD = D(1) - D(2);
    sigma0 = (1 - tanh(tau*(-x_max/2)/x_max))*deltaD/4 + D(2)/2;
    sigmaPhi = (1 - tanh(tau*(x_max/2)/x_max))*deltaD/4 + D(2)/2;
    sigmaKorr = ((2*(sigma0-sigmaPhi)+D(2)-D(1))/(2*x_max))*x - sigma0 + D(1)/2;
    sigma = (1 - tanh(tau*(x-x_max/2)/x_max))*deltaD/4 + D(2)/2 + sigmaKorr;
elseif mode == 2
    % Sigmoidalfunktion nach Lee et al.
    if x/x_max < tau
        omega = 1;
    else
        omega = 1-(x-tau)/(x_max-x_max*tau);
    end
    sigma = ((D(1)-D(2))/2)*(-2*omega^3+3*omega^2) + D(2)/2;
elseif mode == 3
    % Modifizierte Sigmoidalfunktion nach Lee et al.
    if D(1) > D(2)
        omega = 1-(x/x_max)^tau;
    else
        % Dieser Fall wird in der Ausarbeitung nicht berücksichtigt, macht
        % die Funktion aber allgemeingültig. Sie kann damit auch auf den
        % Fall D(2) > D(1) angewendet werden. Dieser kommt praktisch aber
        % nciht vor.
        omega = (1-x/x_max)^tau;
    end
    sigma = ((D(1)-D(2))/2)*(-2*omega^3+3*omega^2) + D(2)/2;
elseif mode == 4  % Im Hauptmenü nicht wählbar/enthalten.
    % Lineare Übergangsfunktion
    sigma = (D(2)-D(1))/(2*x_max)*x+D(1)/2;
end

end
