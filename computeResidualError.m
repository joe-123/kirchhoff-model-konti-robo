function [res] = computeResidualError(u_ini)
    % werte laden
    [tcr,alpha1,alpha2,beta1,beta2,K] = vars();
    % Initialwerte
    u1z = u_ini(3);
    u2z = u_ini(4);
    r_ini = [0 0 0];
    R_ini = [cos(alpha1+beta1*u1z),   -sin(alpha1+beta1*u1z),   0;
             sin(alpha1+beta1*u1z),    cos(alpha1+beta1*u1z),   0;
                        0,                    0,                1];
    theta2_ini = alpha2 + beta2*u2z - (alpha1 + beta1*u1z);
    s_ini = 0;
    %y    = [1 2 3    4 5 6 7 8 9 10 11 12    13 14 15 16    17           18]
    y_ini = [r_ini,   reshape(R_ini,1,9),     u_ini,         theta2_ini,  s_ini];
    y = y_ini;
    s = 0;
    % transition points bestimmen
    p = determineTransitionPoints(tcr, [beta1 beta2]);
    
    % Loop to run consecutive sections of the differential equations.
    for j = 1:length(p)-1
        options2 = odeset('InitialStep',.001);
        u1s=tcr(1,1).tube.u1s;
        u2s=tcr(1,2).tube.u2s;
        %y = y(end,:);
        save('temp','tcr','alpha1','alpha2','beta1','beta2','u1s','u2s','p','K','j');
        [a,b] = ode45(@(s,y) computeDerivative(s,y),[p(1,j) p(1,j+1)],y_ini, options2);
        % ode-funktion läuft in schritten von anfang bis ende des segments
        % ode-funktion sollte array mit Werten zurückgeben, welche die für die KOS / die einzelnen Schritte stehen?!
        % zusätzlich bekommt man auch noch das ergebnis der computeDerivative-Funktion?
        % a = s-wert ; b = y aber es wird dy berechnet und nicht y?!

        % Punkte f�r Abschnitt werden hinuzgef�gt.
        s = [s;a];
        y = [y;b];
        
        % transition conditions
        if j ~= length(p)-1
            % Rotation und Translation soll gleich bleiben, die Steifigkeitsmatrix ändert sich aber
            % --> u1xy, u1z und u2z neu berechnen
            
            % Die Werte mit "_next" im Namen stehen für das nächste Segment.
            % Die Werte ohne "_next" stehen für das aktuelle Segment.
            
            R1 = [y(end,4:6)', y(end,7:9)', y(end,10:12)'];
            theta2 = y(end, 17);
            theta2_dot = y(end,16) - y(end,15);
            R_theta2 = [cos(theta2), -sin(theta2),   0;
                        sin(theta2),  cos(theta2),   0;
                        0,                0,         1];
            R2 = R1*R_theta2;
            K1 = K{1};
            K2 = K{2};
            u1 = y(end, 13:15)';
            u2 = [y(end, 13:14), y(end,16)]';
            u1_star = u1s;
            u2_star = u2s;
            
            K1_next = K1;
            K2_next = K1;
            R1_next = R1;
            R2_next = R2;
            u1_star_next = u1s;
            u2_star_next = u2s;
            
            % Die obenstehenden Werte gelten für den gekrümmten Zustand.
            % Fuer den Fall, dass die tubes im aktuellen/naechsten Segment gerade(1) oder...
            % nicht vorhanden(0) sind, werden die Werte im nachfolgenden korrigiert.
            
            % tube1
            if p(2, j) == 1  % gerade
                u1_star = zeros(3,1);
            end
            
            % tube2
            if p(3, j) == 0  % nicht vorhanden
                u2_star = zeros(3,1);
                K2 = zeros(3,3);
                R2 = zeros(3,3);
                u2 = zeros(3,1);
            elseif p(3, j) == 1  % gerade
                u2_star = zeros(3,1);
            end
            
            %tube1 next
            if p(2, j+1) == 0
                u1_star_next = zeros(3,1);
                K1_next = zeros(3,3);
                R1_next = zeros(3,3);
            elseif p(2, j+1) == 1
                u1_star_next = zeros(3,1);
            end
            
            %tube2 next
            if p(3, j+1) == 0
                u2_star_next = zeros(3,1);
                K2_next = zeros(3,3);
                R2_next = zeros(3,3);
            elseif p(3, j+1) == 1
                u2_star_next = zeros(3,1);
            end

            % Gl.(23) führt auf
            u1_next = inv(R1_next*K1_next + R2_next*K2_next*(R_theta2')) * ...
                         (R1*K1*(u1-u1_star) + R2*K2*(u2-u2_star) + ...
                          R1_next*K1_next*u1_star_next - R2_next*K2_next*(theta2_dot*[0;0;1]-u2_star_next));
            u2_next = R_theta2'*u1_next + theta2_dot*[0;0;1];
            
            y_ini = [y(end,1:12), u1_next(1:2)', u1_next(3), u2_next(3), y(end,17:end)];
        end
    end
    % Fehler berechnen
    u1_tip = y(end, (13:15));
    u2z_tip = y(end,16);
    u1_star_tip = tcr(1).tube.u1s;
    u2_star_tip = tcr(2).tube.u2s;
    moment_res = u1_tip - u1_star_tip';
    res = [moment_res(1:2)'; moment_res(3); (u2z_tip - u2_star_tip(3))]; 
    save('temp','s','y')
end
