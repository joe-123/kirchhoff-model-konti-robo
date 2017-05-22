function [dy] = computeDerivative(s,y)
% hier wird dy mithilfe der kirchhoff rod gleichungen berechnet
    load('temp')
    if p(3, j+1) ~= 0
        K_p = K{1}+K{2};
    else
        K_p = K{1};
    end
    y = y';

    theta2 = y(17);
    
    R_theta1 = [1 0 0;
                0 1 0 ;
                0 0 1];
    R_theta2 = [cos(theta2) -sin(theta2) 0;
                sin(theta2)  cos(theta2) 0;
                0,                0,       1];
    R_theta2_dot = [-sin(theta2) -cos(theta2)  0;
                     cos(theta2) -sin(theta2)  0;
                     0            0            0];
    R = reshape(y(end, 4:12),3,3);
    theta2_dot = y(16) - y(15);
    u1 = y( 13:15);
    u2 = R_theta2'*u1' + theta2_dot*[0;0;1];

    u1_hat = [ 0      -u1(3)    u1(2);
               u1(3)    0      -u1(1);
              -u1(2)    u1(1)    0];
    u2_hat = [ 0      -u2(3)    u2(2);
               u2(3)    0      -u2(1);
              -u2(2)    u2(1)    0];

    % initial shape of tube1 in this segment
    if p(2, j+1) ~= 2
        u1_star = zeros(3,1);
    else
        u1_star = u1s; % gekruemmt
    end
    % initial shape of tube2 in this segment
    if p(3, j+1) ~= 2
        u2_star = zeros(3,1);
    else
        u2_star = u2s;% gekruemmt
    end
    
    % Gl. (5)
    r_dot = R*[0; 0; 1];
    R1_dot = R*u1_hat;

    % Gl. (19) (die Summe wurde der Übersichtlichkeit halber in temp1 und temp2 zerlegt)
    temp1 = R_theta1 *  (u1_hat*K{1} )*(u1'-u1_star);
    if p(3, j+1) ~= 0
        % wird nur berechnet wenn tube2 vorhanden ist:
        temp2 = R_theta2 * (K{2}* (R_theta2_dot*u1') + (u2_hat*K{2} )*(u2-u2_star));
    else
        temp2= zeros(3,1);
    end
    u1xy_dot = -inv(K_p) * (temp1+temp2);
  
    % Gl. (21)
    u1z_dot = ((tcr(1,1).tube.E*tcr(1,1).tube.I) / (tcr(1,1).tube.G*tcr(1,1).tube.J)) * (u1(1)*u1_star(2) - u1(2)*u1_star(1));
    u2z_dot = ((tcr(1,2).tube.E*tcr(1,2).tube.I) / (tcr(1,2).tube.G*tcr(1,2).tube.J)) * (u2(1)*u2_star(2) - u2(2)*u2_star(1));
    if p(3, j+1) == 0
        u2z_dot = y(16);
    end
    
    % Result
    s_dot = y(18);
    dy = [r_dot; reshape(R1_dot,9,1); u1xy_dot(1:2); u1z_dot; u2z_dot; theta2_dot; s_dot];
end
