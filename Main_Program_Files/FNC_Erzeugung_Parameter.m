function [D,L,alpha,beta] = FNC_Erzeugung_Parameter(D0,lim_nTE,Param_Alpha,Param_D,Param_L,mode)

%% Modus: Horsfield et. al. (1981)
nTE = randi([lim_nTE(1),lim_nTE(2)]);
D = zeros(nTE+1,1);
L = D;
alpha = D;
x = Param_D(2);
D_N = (1/nTE)^(1/x)*D0;
if Param_Alpha(3) == 2
    a = (Param_Alpha(2)-Param_Alpha(1)+(Param_Alpha(1)*D_N)/D0)/(D_N^x-D_N*D0^(x-1));
    b = -(Param_Alpha(1)+a*D0^x)/D0;
    fun_alpha = @(d) a*d^x + b*d + Param_Alpha(1);
elseif Param_Alpha(3) == 1
    fun_alpha = @(d) (Param_Alpha(1)/90)*((D0^Param_Alpha(2)-d^Param_Alpha(2))/(D0^Param_Alpha(2)))^(1/Param_Alpha(2));
end
for i = 2 : nTE
    var = Param_D(1);
    while D(i) <= 0 || sum(D.^x) >= D0^x || real(D(i)) ~= D(i)
        % Erzeugung eines normalvertielten Wertes mit dem Erwartungswert
        % D0/nTE 
        D(i) = var*rand + (1/nTE)^(1/x)*D0;
        % Für jeden weiteren Versuch wird die Varianz verringert
        var = var/2;
    end
    % Winkelbetrag abhängig vom Durchmesser
    alpha(i) = fun_alpha(D(i));
    % Elementlänge abhängig vom Durchmesser
    var = Param_L(1);
    if Param_L(3) < Param_L(4)
        while L(i) < Param_L(3)*D(i) || L(i) > Param_L(4)*D(i)
            L(i) = var*randn + Param_L(2)*D(i);
            % Für jeden weiteren Versuch wird die Varianz halbiert
            var = var/2;
        end
    else
        L(i) = Param_L(3)*D(i);
    end
end

% Der letzte Durchmesser wird entsprechend der Fomel nach Horsfield so
% bestimmt, dass der Querschnitt dieser Regel folgt: D0^x = sum(D(i)^x)
D = [D(2:nTE);0];
L = [L(2:nTE);0];
alpha = [alpha(2:nTE);0];
D(nTE) = (D0^x-sum(D.^x))^(1/x);
var = Param_L(1);
if Param_L(3) < Param_L(4)
    while L(nTE) < Param_L(3)*D(nTE) || L(nTE) > Param_L(4)*D(nTE)
        L(nTE) = var*randn + Param_L(2)*D(nTE);
        % Für jeden weiteren Versuch wird die Varianz halbiert
        var = var/2;
    end
else
    L(nTE) = Param_L(3)*D(nTE);
end
% Letzten Winkel berechnen
alpha(nTE) = fun_alpha(D(nTE));
% alpha(nTE) = (2*pi/360)*((abs(Param_Alpha(1))+abs(Param_Alpha(2)))/2)*(1 - D(nTE)/D0);

% Richtung der Winkel bestimmen
for i = 1 : nTE
    if mod(i,2) == 0
        alpha(i) = -alpha(i);
    end
end

% Sortieren der Einträge
DLA = sortrows([D,L,alpha],3,'descend');
D = DLA(:,1);
L = DLA(:,2);
alpha = DLA(:,3);
% beta bestimmen
if mode == 3
    beta = pi*rand - pi/2;
else
    beta = 0;
end