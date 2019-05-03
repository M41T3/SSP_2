%% V2b
clearvars;
clc;

tau = 1;

fg = 1./(2*pi * tau);
f0 = fg * [0.1, 0.5, 1, 2, 10];

H = 1./(1+1i*2*pi*f0*tau);
A =  abs(H);
phi = angle(H);

for f = 1:numel(f0)
   disp("f_0 =  " + num2str(f0(f)) + " - A(" + num2str(f) + ") :" + num2str(A(f)) + " phi("+ num2str(f) +") :" + num2str(phi(f))); 
end

%% V2c
clearvars;
clc;

U = 1;
tau = 1;
T1 = tau * [0.1, 0.5, 1, 2, 10];

for n = 1:numel(T1)
    t_lade = linspace(0, T1(n), 100);
    t_ent = linspace(T1(n), 2*T1(numel(T1)), 100);
    
    uylade = U * (1 - exp(-(t_lade)/tau)); % 0 <= t <= T1
    uyent = U * (1 - exp(-(T1(n))/tau)) * exp((-(t_ent)+T1(n))/tau); % T1 < t <= 2T1

    plots(n) = plot([t_lade, t_ent], [uylade, uyent], 'DisplayName', "T1 = " + num2str(T1(n)));
    hold on;
  
    
end 

hold off;
grid on; grid minor;
title("Systemantwort RC-TP auf rect");
xlabel("t");
ylabel("u(t)");
legend(plots)

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
   figure(1)                                                     %Öffnen eines Plot-Fensters
   plot(t,uy, t,h)                                               %Ausgabe aller Werte               
   xlabel('tau');                                                %Beschriftung x-Achse
   ylabel('U(y) in V');                                          %Beschriftung y-Achse
   title(sprintf("Aufgabe A1: h(t) und Uy(t) pulsbreite tau/%d ", 2^(n-1))); %Titel der Grafik
   grid on;  
   legend('Uy(t) aus P1','h(t)*U*dirac(t)*Tx');
   pause(1.5);                                                   %Pause zum "Halten" des jeweiligen Plots
end


%% A2

clc;
clearvars;


tau = 1;
T1 = tau;
U = 1;
t = linspace(0, 6*tau, 6*256);                                   %Zeitvektor mit 256 Schritten pro tau

ux = U * rectangularPulse((t-T1/2)/T1);
h = (1/tau) * exp(-t/tau);

uy = conv(ux, h, 'same');

plot(t, ux, t, uy, t, h);


