function P = FNC_Drehung3xz0(P,n,psi)

n = n/norm(n);

%% Rotationsmatrix
mc = 1 - cos(psi);
R  = [cos(psi) + n(1) * n(1) * mc,    n(1) * n(2) * mc,                n(2) * sin(psi); ...
      n(1) * n(2) * mc,               cos(psi) + n(2) * n(2) * mc,     -n(1) * sin(psi); ...
      -n(2) * sin(psi),               n(1) * sin(psi),                 cos(psi)];

%% Rotation
P = P*R;
    
end

