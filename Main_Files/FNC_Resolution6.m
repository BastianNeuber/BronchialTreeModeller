function [N] = FNC_Resolution6(n)
% Die Funktion bestimmt aus einem Eingabewert die nächsthöhere ganze gerade
% Zahl, die nicht durch 4 Teilbar und größer oder gleich 6 ist.

n = ceil(n);
if mod(n,4) == 0
    N = n + 2;
elseif mod(n,4) ~= 2
    N = n + mod(n,4);
else
    N = n;
end

end