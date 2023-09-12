function [estimanted_pos]=estimated_coordinate_RSS_D2D(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,Blockage_real,Blockaget_est,nu_real,nu_est,sigma)

    min_allXY = min(X,[],1);
    max_allXY = max(X,[],1);

    %% non-changing parameters
    frequency = 28; %GHz
    R = 1; %reference distance for path loss
    Pt = db2pow(10)/1000; % Transmit Power 10 db in [W]
    G_rx_lin = db2pow(5.57); % 5.57 dBi Receive gain in linear scale
    G_lin = db2pow(14.58); % 14.58 dBi Transmit gain in linear scale
  
    %% noise parameters (changing)
    mu = 0;
%     sigma = 4; % 4dB for LOS, 7.82 NLOS
    w_dB_real = normrnd(mu,sigma); % noise for NLOS
    w_dB_est = normrnd(mu,sigma);

    %% Noise power and interference
    B = 10e7; %bandwidth in Hz 100 MHz
    noisefiguredB = 10;
    % compute the noise power in dB
    sigma2dBm = -174 + 10*log10(B) + noisefiguredB;
    sigma2 = db2pow(sigma2dBm);
    Is = 10*sigma2dBm*Pt;
%      Is=0;

    %% 1 meter PL and 1 Meter RSSI in dB
    PL_real = db2pow(20*log10(frequency) + (nu_real*10)*log10(R) + Blockage_real - Is); %PL(d0) in Watts d0=R
    PL_est = db2pow(20*log10(frequency) + (nu_est*10)*log10(R) + Blockaget_est - Is); %PL(d0) in Watts
    
    Measured_power = real(pow2db((Pt*G_lin*G_rx_lin)/PL_real)); % 1 Meter RSSI in dB (watt to dB)    
    Measured_power_est = real(pow2db((Pt*G_lin*G_rx_lin)/PL_est)); % 1 Meter RSSI in db (watt to dB)    
   
    %% measure real RSS  
    d_real(1) = real(sqrt((X_a1 - realX)^2 + (Y_a1 - realY)^2)); %real distance 
    d_real(2) = real(sqrt((X_a2 - realX)^2 + (Y_a2 - realY)^2)); %real distance 
    d_real(3) = real(sqrt((X_a3 - realX)^2 + (Y_a3 - realY)^2)); %real distance 
    d_real(4) = real(sqrt((X_a4 - realX)^2 + (Y_a4 - realY)^2)); %real distance 
    
    RSS = round(Measured_power - 10.*nu_real.*log10(d_real) + w_dB_real); %we round RSS to the interger value
    
    %% distance estimation
    %RSS is real, others are unknown
    result = (Measured_power_est - RSS + w_dB_est)/(10*nu_est); 
    
    d = estimate_distance(result);
    
    estimanted_pos = real(trilateration_estimation(X, d, realX, realY, min_allXY(1), max_allXY(1), min_allXY(2), max_allXY(2)));

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

