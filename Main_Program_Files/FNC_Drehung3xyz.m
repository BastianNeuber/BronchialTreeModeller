function P = FNC_Drehung3xyz(A, B, P, alpha)
%% Die Funktion stammt aus dem Mathworks FileExchange: https://de.mathworks.com/matlabcentral/fileexchange/66446-rotation-matrix?focused=8946659&tab=function

%% Verschiebung in den Ursprung des Koordinatensystems
P = P - ones(size(P,1),1)*B;

%% Rotationsachse
n = A - B;
if norm(n) ~= 0
    n = n/norm(n);
    s = sin(alpha); 
    c = cos(alpha); 

    %% 3D rotation matrix:
    x  = n(1);
    y  = n(2);
    z  = n(3);
    mc = 1 - c;
    R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
        x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
        x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];

    %% Rotation
    P = P*R;
end

%% RÜckverschiebung
P = P + ones(size(P,1),1)*B;

end