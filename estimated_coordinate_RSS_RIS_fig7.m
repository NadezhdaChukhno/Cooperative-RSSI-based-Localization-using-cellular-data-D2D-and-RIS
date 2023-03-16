function [estimanted_pos,snr]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN)

    min_allXY=min(X,[],1);
    max_allXY=max(X,[],1);

    %% non-changing parameters
%     frequency=28; %GHz
    frequency=200; %GHz
    lambda =0.01071; % wavelength 
    %L = fspl(R,lambda) returns the free space path loss in decibels for a waveform with wavelength lambda propagated over
    %a distance of R meters. The minimum value of L is zero, indicating no path loss.
    R=1; %reference distance
    Pt_dBm = 20; % Transmit Power in dBm
    Pt = 10^(Pt_dBm/10)/1000; % Transmit Power  [W]
    G_rx=5.57; %dBi Receive gain 
    G_rx_lin=10^(G_rx/10); % Reveive gain in linear scale
    G=14.58; %dBi
    G_lin=10^(G/10); %dBi

    %% changing parameters 
    Blockage_real=32.4; % blockage real
    Blockaget_est=32.4; % blockage estimated
    nu_real=2.1; %path loss expanent; %LOS and NLOS path loss exponent
    % gaussian PLE
    mu_pl=2;
    sigma_pl=3;
%     nu_est=abs(normrnd(mu_pl,sigma_pl));
      nu_est=2.1; 

    %% noise parameters (changing)
    % Gaussian USED for NLOS
    mu=0;
    sigma=4; % 4dB for LOS, 7,82 NLOS
    w_dB_real = normrnd(mu,sigma); % noise for NLOS
    w_dB_est = normrnd(mu,sigma);

    % left skewed USED for LOS
    skew=2;
    kurt=6;
    %w_dB_real  = pearsrnd(mu,sigma,skew,kurt);
    N0=10^(-174/10); %Power spectral density of noise, N0 

    W_snr=1000000000; %1GHz (Hz)
    
    
    Gamma=1; %reflection gain from RIS
    %% 1 meter PL and 1 Meter RSSI in dB
%     PL_db_real = 20*log10(frequency)+(nu_real*10)*log10(R)+Blockage_real; %PL(d0) in dB
%     PL_real=10^(PL_db_real/10); %PL(d0) in Watts
    
    PL_db_real_RD = 20*log10(frequency)+(nu_real*10)*log10(R)+Blockage_real; %PL(d0) in dB
    PL_real_RD=10^(PL_db_real_RD/10); %PL(d0) in Watts
    
    PL_db_real_SR = 20*log10(frequency)+(nu_real*10)*log10(R)+Blockage_real; %PL(d0) in dB
    PL_real_SR=10^(PL_db_real_RD/10); %PL(d0) in Watts
    
    PL_real=((sqrt((1/(PL_real_RD*PL_real_SR))))*NN)^(-2);
    
%     PL_db_est_RD = 20*log10(frequency)+(nu_est*10)*log10(R)+Blockaget_est; %PL(d0) in dB
%     PL_est_RD=10^(PL_db_est_RD/10); %PL(d0) in Watts
    PL_est_RD=10^(2*log10(frequency)+Blockaget_est*0.1)*R^(nu_est);
    
%     PL_db_est_SR = 20*log10(frequency)+(nu_est*10)*log10(R)+Blockaget_est; %PL(d0) in dB
%     PL_est_SR=10^(PL_db_est_SR/10); %PL(d0) in Watts
    PL_est_SR=10^(2*log10(frequency)+Blockaget_est*0.1)*R^(nu_est);

    PL_est=((sqrt((1/(PL_est_RD*PL_est_SR))))*NN)^(-2);
    
    Measured_power_W=(Pt*G_lin*G_rx_lin)/PL_real; % 1 Meter RSSI in Watt
    Measured_power=real(10*log10(Measured_power_W)); % 1 Meter RSSI in dB
        
    Measured_power_W_est=(Pt*G_lin*G_rx_lin)/PL_est; % 1 Meter RSSI in Watt
    Measured_power_est=real(10*log10(Measured_power_W_est)); % 1 Meter RSSI in dB

    SNR_watt=((Pt*G_rx_lin*G_lin*Gamma)/(N0*W_snr*PL_real)); % I disagree with formula (15) from the file
    snr=real(10*log10(1000*SNR_watt)); 
        
    
    
    %% measure real RSS  

    d_real_RD(1)=real(sqrt(( X_a1- realX)^2 + ( Y_a1- realY)^2 )); %real distance 
    d_real_RD(2)=real(sqrt(( X_a2- realX)^2 + ( Y_a2- realY)^2 )); %real distance 
    d_real_RD(3)=real(sqrt(( X_a3- realX)^2 + ( Y_a3- realY)^2 )); %real distance 
    d_real_RD(4)=real(sqrt(( X_a4- realX)^2 + ( Y_a4- realY)^2 )); %real distance 
    
    d_real_SR(1)=real(sqrt(( X_a1_BS- X_a1)^2 + ( Y_a1_BS- Y_a1)^2 )); %real distance 
    d_real_SR(2)=real(sqrt(( X_a2_BS- X_a2)^2 + ( Y_a2_BS- Y_a2)^2 )); %real distance 
    d_real_SR(3)=real(sqrt(( X_a3_BS- X_a3)^2 + ( Y_a3_BS- Y_a3)^2 )); %real distance 
    d_real_SR(4)=real(sqrt(( X_a4_BS- X_a4)^2 + ( Y_a4_BS- Y_a4)^2 )); %real distance 
    
    RSS= round(Measured_power-10.*nu_real.*log10(d_real_RD+d_real_SR)+w_dB_real);
%     RSS= round(Measured_power-10.*nu_real.*log10(d_real_RD).*log10(d_real_SR)+w_dB_real);
    
    %% distance estimation
    %RSS is real, others are unknown
    result=(Measured_power_est-RSS+w_dB_est)/(10*nu_est); 
    
    d_ALL=estimate_distance(result);
    
    d=d_ALL-d_real_SR';
    
    estimanted_pos=real(trilateration_estimation(X, d, realX, realY, min_allXY(1), max_allXY(1), min_allXY(2), max_allXY(2)));
    %  X - Matrix containing the coordinates for each AP position
    %  d - Estimated distances to each AP, respectively
    %  realX - x coordinate of the known position
    %  realY - y coordinate of the known position
    %  min_X - minimum X coordinate of the dataset
    %  max_X - maximum X coordinate of the dataset
    %  min_Y - minimum Y coordinate of the dataset
    %  max_Y - maximum Y coordinate of the dataset
    %  OUTPUTS: Estimated position (x, y)
    

end

