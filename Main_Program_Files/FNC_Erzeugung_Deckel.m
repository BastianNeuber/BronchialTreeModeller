function [P] = FNC_Erzeugung_Deckel(n,Di,Da)

P = [];
for i = ceil((Di/Da)*n) : n
    for j = 1 : n
        P = [P;[(Da/2)*(i/n)*cos(j/n*2*pi),(Da/2)*(i/n)*sin(j/n*2*pi),0]];
    end
end

end

