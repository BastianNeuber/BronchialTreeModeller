clc
Param = app.UITable.Data;
n = app.PunktabstandEditField.Value/1000;
P = [];
Pa = [];
Ps = [];
lastrow = 0;

%% Innenform konstruieren

for i = 1 : size(Param,1)
    if isequal(Param(i,2),{'Kegel'})
        D = str2double(Param(i,6))/2;
        alpha = deg2rad(str2double(Param(i,5)));
        L = D / sin(alpha);
        N = APPFNC_ceil2(L/n);  % N(i) = APPFNC_ceil2(L/n);
        p = APPFNC_Kegel(N(i),L,alpha);
        lastrow = lastrow + size(p,1);
    elseif isequal(Param(i,2),{'Kreisbogen'})
        r = str2double(Param(i,3));  % 1. Radius des Kreises
        phi = str2double(Param(i,5))*pi/180;  % Winkel, den der Kreisbogen umschließt
        S = P(lastrow,:);
        p = APPFNC_Kreisbogen(n,S,r,alpha,phi);  % [p,N(i)] = APPFNC_Kreisbogen(n,S,r,alpha,phi);
        alpha = alpha - phi*sign(r);
        lastrow = lastrow + size(p,1);
    elseif isequal(Param(i,2),{'Strecke'})
        L = str2double(Param(i,4));  % Länge der Strecke
        beta = deg2rad(str2double(Param(i,5)));  % Winkel, den die Strecke mit der z-Achse bildet
        S = P(lastrow,:);
        p = APPFNC_Strecke(n,S,L,beta);  % [p,N(i)] = APPFNC_Strecke(n,S,L,beta);
        lastrow = lastrow + size(p,1);
        alpha = beta;
    end
    P = [P;p];
end


%% Außenlinie konstruieren

if app.Wandstaerke_CB.Value == true
    Pa = zeros(size(P));
    WS = app.Wandstaerke_EF.Value;

    % Horizont
    h = min(ceil(WS/n),floor(size(P,1)/2));
    
    % Schleife über alle Punkte, bei denen nur der Horizont vorwärts gültig ist
    % (1 bis i+h)
    for i = 1 : h  % min(2*h,size(P,1)-h)
        if P(i+h,2)-P(1,2) > 0
            lambda = atan((P(i+h,1)-P(1,1))/(P(i+h,2)-P(1,2)));
        elseif P(i+h,2)-P(1,2) < 0
            lambda = atan((P(i+h,1)-P(1,1))/(P(i+h,2)-P(1,2))) + pi/2;
        else
            lambda = pi/2;
        end
        Pa(i,:) = P(i,:) + [sin(lambda+pi/2),cos(lambda+pi/2)]*WS;  % *sign(90.000001+str2double(Param{j,5}));
    end

    % Schleife über alle Punkte, bei denen der Horizont in beide Richtungen
    % gültig ist (i-h bis i+h)
    for i = h+1 : size(P,1)-h
        if P(i+h,2)-P(i-h,2) > 0
            lambda = atan((P(i+h,1)-P(i-h,1))/(P(i+h,2)-P(i-h,2)));
        elseif P(i+h,2)-P(i-h,2) < 0
            lambda = atan((P(i+h,1)-P(i-h,1))/(P(i+h,2)-P(i-h,2))) + pi/2;
        else
            lambda = pi/2;
        end
        Pa(i,:) = P(i,:) + [sin(lambda+pi/2),cos(lambda+pi/2)]*WS;  % *sign(90.000001+str2double(Param{j,5}));
    end

    % Schleife über alle Punkte, bei denen der Horizont nur rückwärts gültig
    % ist (i-h bis Ende)
    I = size(P,1);
    for i = size(P,1)-h+1 : size(P,1)
        if P(I,2)-P(i-h,2) > 0
            lambda = atan((P(I,1)-P(i-h,1))/(P(I,2)-P(i-h,2)));
        elseif P(I,2)-P(i-h,2) < 0
            lambda = atan((P(I,1)-P(i-h,1))/(P(I,2)-P(i-h,2))) + pi/2;
        else
            lambda = pi/2;
        end
        Pa(i,:) = P(i,:) + [sin(lambda+pi/2),cos(lambda+pi/2)]*WS;  % *sign(90.000001+str2double(Param{j,5}));
    end
else
    Pa = [];
end


%% Außen- und Innenwand verbinden

if app.EndstueckSchliessen_CB.Value == true
    N_es = ceil(norm(P(end,:)-Pa(end,:))/n);
    Ps = zeros(N_es-1,2);
    for i = 1 : N_es-1
        Ps(i,:) = P(end,:) + (Pa(end,:)-P(end,:))*i/N_es;
    end
else
    Ps = [];
end


%% Punktewolke speichern und plotten

save('Settings\Punkte_Endstueck.mat','P','Pa','Ps');

% Plot update
if size(app.UITable.Data,1) ~= 0
    scatter(app.UIAxes,P(:,1),P(:,2),8,'filled');
    if app.Wandstaerke_CB.Value == true
        hold(app.UIAxes);
        scatter(app.UIAxes,Pa(:,1),Pa(:,2),8,'filled');
        hold(app.UIAxes);
    end
    if app.EndstueckSchliessen_CB.Value == true
        hold(app.UIAxes);
        scatter(app.UIAxes,Ps(:,1),Ps(:,2),8,'filled');
        hold(app.UIAxes);
    end
    axis(app.UIAxes,'equal')
else
    scatter(app.UIAxes,[],[]);
end


%% Anzeige einstellen
% P = [P;Pa;Ps];
% % % Ticklänge
% % l_Tick = round(max(P(:,1))/8,2,'significant');
% % x_Tick = 0:l_Tick:max(P(:,1));
% % XTick(app.UIAxes,x_Tick);
% if max(P(:,2))*2 < max(P(:,1))
%     xlim(app.UIAxes,[0,1.1*max(P(:,1))])
%     ylim(app.UIAxes,[min(P(:,2)),min(P(:,2))+2*max(P(:,1))])
% else
%     
% end
% % axis(app.UIAxes,'equal')