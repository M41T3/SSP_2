%% A1

clc;
clearvars;

% Verlauf der Sprungantwort bei Ansteuerung mit Approximaxion von einem
% Dirac

tau = 1;
U = 1;
Tx = 1;

t = linspace(0, 6*tau, 6*256);                                   %Zeitvektor mit 256 Schritten pro tau

% Berchnen des Ausgangssignals mit Dirac mit variabler Pulslänge
for n = 1:8
    
    Tn = tau / (2^(n-1));                                        %Verkleinerung der Impulsbreite
    
    %Uy aus Praktikum 1
    Utn = U * Tx * (1/Tn) * (1-exp(-Tn/ tau) );                  %Spannung bei Uy(tn)
    
    uy = (t<=Tn) .* U .* Tx .* (1/Tn) .* (1-exp(-t/ tau) )...    %Ladefunktion
        + (t > Tn).* Utn .* exp(-(t - Tn) / tau);                %Entladefunktion <- richtig baba!
    
    h = (1/tau) .* exp(-t/tau);
    %h = heaviside(t) .* (1/tau) .* exp(-t/tau);
    
   %Grafikausgabe
   figure(1)                                                     
   plot(t,uy, t,h)                                               %Ausgabe aller Werte               
   xlabel('tau');                                                
   ylabel('U(y) in V');                                          
   title(sprintf("Aufgabe A1: h(t) und Uy(t) pulsbreite tau/%d ", 2^(n-1)));
   grid on;  
   legend('Uy(t) aus P1','h(t)*U*dirac(t)*Tx');
   pause(1.5);                                                   
end


%% A2

clc;
clearvars;

% Faltung von ux * h(t)
tau = 1;
T1 = tau;
U = 1;


% Faltung
T = (6*tau)/(1536);                                              % 1536 Messwerte pro tau
k = 1:1536;

tux = (k-1).*T;                                                  % Zeitschritte für conv


ux = heaviside(tux) .* U - heaviside(tux-T1) .* U;               % heaviside, da rect(0) = 1/2
h = heaviside(tux) .* (1/tau) .* exp(-(tux)/tau);    
uy = conv(ux, h) .*T;                                           % Näherungsfunktion Faltung

% Aus V1
t = linspace(0, 6*tau, 256*6);                                   %Zeitvektor mit 256 Schritten pro tau
uyt = (t <= T1) .* U .* (1 - exp(-t/tau)) ...
    + (t > T1) .* U .* (1 - exp(-T1/tau)) .* exp(-(t-T1)/tau);

nT = 1:(2*1536)-1;


figure(1)                                                        %Öffnen eines Plot-Fensters
plot(nT.*T, uy, t, uyt);
xlabel('T1 = tau');                                              %Beschriftung x-Achse
ylabel('Uy in V');                                               %Beschriftung y-Achse
title(sprintf("Aufgabe A2: f = Uy(t) "));                        %Titel der Grafik
grid on; grid minor; 
legend('Uy(t) aus Faltung','Uy(t) Th. Erg. aus V1');


%% A3

clc;
clearvars;

tau = 1;
fg = 1/(2*pi*tau);
f0 = fg .*[1/10, 1/2, 1, 2, 10];
U = 1;
res = 2000;
for i = 1:1
    T = 10*tau + 2/f0(i);
    Tconv = T/(res);                                             % Abstand Zeitschritte für Faltung

    k = 1:res;
    tux = (k-1).*Tconv;                                          % Zeitschritte für conv

    t = linspace(0, 2*T, res);
    ux = (heaviside(t) .* U - heaviside(t-T) .* U) .* sin(2.*pi.*f0(i).*t);

    uxf = (heaviside(tux) .* U - heaviside(tux-T) .* U) .* sin(2.*pi.*f0(i).*tux);
    h = heaviside(tux) .* (1/(tau)) .* exp(-(tux)/(tau));    
    uy = conv(uxf, h) .* Tconv;                                  % Näherungsfunktion Faltung

    nT = 1:(2*res)-1;



    figure(1);                                                   %Öffnen eines Plot-Fensters
    plot(t, ux ,nT.*Tconv, uy);          
    xlabel('T1 = tau');                                          %Beschriftung x-Achse
    ylabel('Uy in V');                                           %Beschriftung y-Achse
    title(sprintf("Aufgabe A3: f = Uy(t) mit f0=%.2d", f0(i)));  %Titel der Grafik
    grid on; grid minor;
    legend('Ux(t)','Uy(t)');
    pause(1.5);
end
