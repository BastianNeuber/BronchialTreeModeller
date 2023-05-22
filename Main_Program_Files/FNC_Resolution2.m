function [N] = FNC_Resolution2(n)
% Die Funktion bestimmt aus einem Eingabewert die nächsthöhere ganze gerade
% Zahl, die nicht durch 4 Teilbar und größer oder gleich 6 ist.

n = ceil(n);
if n <= 0
    N = 2;
elseif mod(n,2) == 1
    N = n + 1;
else
    N = n;
end

end