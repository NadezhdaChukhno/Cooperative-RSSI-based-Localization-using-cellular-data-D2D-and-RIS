%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig.2 SNR comparison of RIS and relays                     %
% Article: [name]                                                                %
% Download article: [link]                                                       %
% This is version 1.0 (Last edited: 2022-07-13)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
%% parameters 
frequency=28; % frequency [GHz]
N0=10^(-174/10); % power spectral density of noise, N0 
Pt_dBm = 23.010299957; % transmit power [dB]
Pt = 10^(Pt_dBm/10)/1000; % transmit power[W]
W=1;% bandwidth [GHz]
W_snr=1000000000; % bandwidth [Hz]
G_rx=5.57; % received gain [dBi]
G_rx_lin=10^(G_rx/10); % received gain [linear scale]
G =14.58; % transmit gain [dBi] (32 antenna elements)
G_lin =10^(G/10);  % transmit gain in linear scale
distance=100; % transmission distance [m] between tranmitter and reciever    

%% RIS parameters
Gamma=1; % reflection gain from RIS
NN = [256,1024,2048]; % number of reflective elements 
a=0.1:0.1:0.9; %at which distance from BS/UE we need to place RIS (0y axis on Figure 2)

%% simulation
for i=1:length(a)
    for n=1:length(NN)
        total_loss_nB_SR(i)=10^(2*log10(frequency)+3.24)*(a(i)*distance)^(2.1); % linear non blocked path loss 
        total_loss_nB_RD(i)=10^(2*log10(frequency)+3.24)*((1-a(i))*distance)^(2.1); % linear non blocked path loss
        total_loss_nB(i)=   10^(2*log10(frequency)+3.24)*(distance)^(3.19); %NLOS no bl
        total_loss_nB2(i)=  10^(2*log10(frequency)+4.74)*(distance)^(2.1);% LOS bl
        
        % total path loss
        LRIS(i,n)=((sqrt((1/(total_loss_nB_SR(i)*total_loss_nB_RD(i)))))*NN(n))^(-2); % formula (23) from article

        % snr
        SNR_watt(i,n)=((Pt*G_rx_lin*G_lin*Gamma)/(N0*W_snr*LRIS(i,n))); 
        snr(i,n)=real(10*log10(1000*SNR_watt(i,n)));   %W is in Hz   
         
        %snr without RIS NLOS no blockage
        SNR_watt_woRIS(i,n)=(Pt*G_rx_lin*G_lin)/(N0*W_snr*total_loss_nB(i));
        snr_woRIS(i,n)=real(10*log10(1000*SNR_watt_woRIS(i,n)));    %W is in Hz   
        
        %snr without RIS LOS blockage
        SNR_watt_woRIS2(i,n)=(Pt*G_rx_lin*G_lin)/(N0*W_snr*total_loss_nB2(i)); 
        snr_woRIS2(i,n)=real(10*log10(1000*SNR_watt_woRIS2(i,n)));    %W is in Hz   
    end
end

%% Plot the curves
        
name_NN1 = strcat('256 RIS',{' '});
name_NN2 = strcat('1024 RIS',{' '});
name_NN3 = strcat('2048 RIS', {' '});
name_W = strcat('Without RIS NLoS no blockage', {' '});  
name_W2 = strcat('Without RIS LoS blockage', {' '});  

 figure(1)
 hold on
 plot(a, snr(:,1), 'DisplayName', name_NN1{1}, ...
    'LineWidth', 1, 'Marker', 'square');

 plot(a, snr(:,2), 'DisplayName', name_NN2{1}, ...
    'LineWidth', 1, 'Marker', '>');

 plot(a, snr(:,3), 'DisplayName', name_NN3{1}, ...
    'LineWidth', 1,   'Marker', 'o');

 plot(a, snr_woRIS(:,1), 'DisplayName', name_W{1}, ...
    'LineWidth', 1,   'Marker', 'o');
 plot(a, snr_woRIS2(:,1), 'DisplayName', name_W2{1}, ...
    'LineWidth', 1,   'Marker', 'o');

legend('show'); grid on
ylabel('SNR, dBm')
xlabel('d_{SR} and d_{RD} distance relationship')