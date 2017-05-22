clc
clear vars
close all
pkg load odepkg             % <<<<<<<<<<<<<<<<<<<<<<<<<<< das ist nur für octave !!!-------------------


% --->> values are in "vars.m" !!!
[tcr,alpha1,alpha2,beta1,beta2,K] = vars();
u_initial = [0 0 0 0];

options1 = optimset('TolFun', 1e-6);
if (beta1 < 0.2 && beta2 < 0.1)
[q,w] = fsolve(@computeResidualError, u_initial, options1);
load('temp'); % zwischengespeicherte werte laden
g = [y(:,4:6),  zeros(size(y,1),1)    ,y(:,7:9),    zeros(size(y,1),1), y(:,10:12), zeros(size(y,1),1), y(:,1:3),   ones(size(y,1),1)];
else
    % tubes komplett eingefahren
    g=[1,zeros(1,4),1,zeros(1,4),1,zeros(1,4),1];
    s=0;
end
drawTCR(g,s,tcr,[-beta1,-beta2])
