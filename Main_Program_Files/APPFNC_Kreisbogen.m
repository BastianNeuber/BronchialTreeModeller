function p = APPFNC_Kreisbogen(n,S,r,alpha,phi)
% n = Punktabstand
% s = Startpunkt
% r = Kreisradius
% alpha = Startwinkel, den der Kreisbogen am Kontaktpunkt hat
% phi = Winkel, der durch den Kreisbogen eingeschlossen wird

% Punkteanzahl berechnen
N = abs(ceil((r*phi/n)));

% Winkelvektor berechnen
M = S + r*[-cos(alpha) , sin(alpha)];
p = zeros(N,2);
Phi = linspace(alpha+sign(r)*pi/2,alpha+sign(r)*(pi/2-phi),N);

% Ellipsenbogen berechnen
for i = 1 : N
    p(i,:) = M + abs(r)*[sin(Phi(i)),cos(Phi(i))];
end

end