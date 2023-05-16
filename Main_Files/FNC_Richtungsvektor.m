function P = FNC_Richtungsvektor(A, B)

%% Rotationsachse
n = [B(2) - A(2); -B(1) + A(1); 0];
if norm(n) ~= 0
    n = n/norm(n);

    %% Rotationswinkel
    if B(3) > A(3)
        psi =         asin(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A));
    elseif B(1) > A(1)
        psi = (pi/2 + acos(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A)))*sign(B(1)-A(1));
    elseif B(1) < A(1)
        psi = (pi   + asin(sqrt((B(1)-A(1))^2 + (B(2) - A(2))^2)/norm(B-A)))*sign(B(1)-A(1));
    end

    %% Rotationsmatrix
    epsilon_r = 1 - cos(psi);
    R  = [cos(psi) + n(1) * n(1) * epsilon_r,    n(1) * n(2) * epsilon_r,                n(2) * sin(psi); ...
          n(1) * n(2) * epsilon_r,               cos(psi) + n(2) * n(2) * epsilon_r,     -n(1) * sin(psi); ...
          -n(2) * sin(psi),                      n(1) * sin(psi),                        cos(psi)];    
    
    %% Rotation
    P = [0,0,1]*R;
else
    P = [0,0,1];
end

end

