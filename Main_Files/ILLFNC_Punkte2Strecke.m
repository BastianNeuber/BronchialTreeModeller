function [P] = FNC_Punkte2Strecke(A,B,n)
%Die Funktien gibt eine Matrix mit n Punkten aus, die auf einer Strecke
%liegen. Die Strecke beginnt in A und endet in B.

P = zeros(n,size(A,2));
for i = 1 : n
    P(i,:) = i/n*(B-A) + A;
end
P = [A;P];

end

