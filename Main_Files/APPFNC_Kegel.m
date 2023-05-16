function [p] = APPFNC_Kegel(N,L,alpha)
p = zeros(N+1,2);
p(1,2) = -L*cos(alpha);
if cos(alpha) > 1e-9
    for i = 2 : N+1
        p(i,:) = [L*((i-1)/N)*sin(alpha),L*(i-1)/N*cos(alpha)-L*cos(alpha)];
    end
else
    for i = 2 : N+1
        p(i,:) = [L*(i-1)/N,0];
    end
end
end