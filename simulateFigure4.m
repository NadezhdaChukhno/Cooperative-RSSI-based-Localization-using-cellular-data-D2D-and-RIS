%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig.4 Visualisation of distance-dependent noise variances, 
% 3GPP UMi Street Canyon LoS/nLoS                   %
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 4.0 (Last edited: 2024-02-20)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ours
%UMi 
d_array = 0:5:600; % min distance:step:max distance
syms x
dLOS10 = solve(0.1 == 18/x + exp(-x/36)*(1 - 18/x), x);
dLOS20 = solve(0.2 == 18/x + exp(-x/36)*(1 - 18/x), x);
dLOS30 = solve(0.3 == 18/x + exp(-x/36)*(1 - 18/x), x);
dLOS70 = solve(0.7 == 18/x + exp(-x/36)*(1 - 18/x), x);

for i = 1:length(d_array)
    d = d_array(i);
    fc_GHz = 28;
%     d = 400;
    FSPL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d); % dB 
    FSPL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d); % dB 
    
    if d <= 18
        P_LOS_UMI = 1;
    else
        P_LOS_UMI = 18/d + exp(-d/36)*(1 - 18/d);
    end
    
    PL_theor(i) = FSPL_LOS * P_LOS_UMI + FSPL_NLOS * (1 - P_LOS_UMI);
    % noise in the system considering bandwidth
    Pn = db2pow(-173.10 + 10*log10(100*10^6)); % noise linerar 354.81 K noise temperature gives -173.10 dBm/Hz
    Pt = db2pow(20)/1000;
    noise_var(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_theor(i)) / ( (100*10^6)^2 *(Pt/Pn) ) );
    
    %[here]
    %% empirical 10%
    % let's take 10% LOS prob. Then disctance is approx 200 in LOS
    if d>dLOS10
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS10);
        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS10);
    else
        PL_NLOS = 0;
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
    end
    PL_emp10(i) = PL_LOS + PL_NLOS;
    noise_var_emp10(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp10(i)) / ( (100*10^6)^2 *(Pt/Pn) ));

    %% empirical 20%
    if d>dLOS20
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS20);
        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS20);
    else
        PL_NLOS = 0;
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
    end
    PL_emp20(i) = PL_LOS + PL_NLOS;
    noise_var_emp20(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp20(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
    %% empirical 30%
    if d>dLOS30
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS30);
        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS30);
    else
        PL_NLOS = 0;
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
    end
    PL_emp30(i) = PL_LOS + PL_NLOS;
    noise_var_emp30(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp30(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
    %% empirical 70%
 
    
    if d>dLOS70
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS70);
        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS70);
    else
        PL_NLOS = 0;
        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
    end
    PL_emp70(i) = PL_LOS + PL_NLOS;
    noise_var_emp70(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp70(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
end

plot(d_array,PL_theor,'LineWidth', 1.5, 'Color', [44/255 160/255 44/255])
hold on
plot(d_array,PL_emp10,'LineWidth', 1.5, 'Color', [0 0 1])
hold on
plot(d_array,PL_emp20,'LineWidth', 1.5, 'Color', [1 0 0])
hold on
plot(d_array,PL_emp30,'LineWidth', 1.5, 'Color', [0.4 0.2 0.8])
hold on
plot(d_array,PL_emp70,'LineWidth', 1.5, 'Color', [0.2 0.6 1])
ylabel('Path loss, dB')
xlabel('Distance between BS and MT, d_{2D}, m')
legend('Theoretical','Empirical, guaranteed 10% LOS','Empirical, guaranteed 20% LOS','Empirical, guaranteed 30% LOS','Empirical, guaranteed 70% LOS')
grid on

figure
plot(d_array,noise_var,'LineWidth', 1.5, 'Color', [44/255 160/255 44/255])
hold on
plot(d_array,noise_var_emp10,'LineWidth', 1.5, 'Color', [0 0 1])
hold on
plot(d_array,noise_var_emp20,'LineWidth', 1.5, 'Color', [1 0 0])
hold on
plot(d_array,noise_var_emp30,'LineWidth', 1.5, 'Color', [0.4 0.2 0.8])
hold on
plot(d_array,noise_var_emp70,'LineWidth', 1.5, 'Color', [0.2 0.6 1])
ylabel('Noise variance, dB')
xlabel('Distance between BS and MT, d_{2D}, m')
legend('Theoretical','Empirical, guaranteed 10% LOS','Empirical, guaranteed 20% LOS','Empirical, guaranteed 30% LOS','Empirical, guaranteed 70% LOS')
grid on





