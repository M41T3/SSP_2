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
