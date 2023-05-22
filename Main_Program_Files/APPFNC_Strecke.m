function p = APPFNC_Strecke(n,S,L,beta)

N = APPFNC_ceil2(L/n);  % Anzahl der Punkte
p = zeros(N,2);  % Punktematrix vorbereiten
for i = 1 : N  % Schleife über alle Punkte
    p(i,:) = S + i/N*[sin(beta),cos(beta)]*L;
end

end