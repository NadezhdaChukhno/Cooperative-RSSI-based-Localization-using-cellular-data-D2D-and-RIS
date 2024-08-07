clear all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig. 6 Theoretical lower bounds, equal placement of RISs and relays                                                     %
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 4.0 (Last edited: 2024-02-20)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use mode 3
%% RSS  location estimation
%% parameters
% prompt_number_of_evaluation_points = "Enter number of evaluation points (M): ";
% number_of_evaluation_points = input(prompt_number_of_evaluation_points);
%number_of_evaluation_points = 30;
% prompt_runs = "Enter number of experiment repetitions (rep): ";
% runs = input(prompt_runs);
% runs = 10;
%% BS/RIS/USER/D2D coordinates
% 4 coordinates of known anchors
number_of_anchors = 4;
Rd = 100; % 100 x 100 square
mode_array = [1,2,3];
% mode = 3;
%1 - cellular    
%2 - D2D
%3 - RIS
NN = 1024; % Number of reflective elements
%% D2D mode switching
prompt_d2d_mode = "Choose D2D user distribution? Options are 1,2,3,4: ";
mode_d2d = input(prompt_d2d_mode);
%mode_d2d = 1;
%1 - random 4 points within the considered area
%2 - random 4 points near real point within the radius of 20 m
d0 = 20; %radius of 20 m
%3 - equa3l placement of RISs and relays
%4 - 4 points on the corners on variable (5,20) distance [not used in the
%paper]
%fixed random position of d2d anchors for mode 1
r2 = Rd*sqrt(rand(number_of_anchors,1));
theta2 = (pi/2)*rand(number_of_anchors,1);
X_BF = (r2.*cos(theta2))';
Y_BF = (r2.*sin(theta2))';  
%% Cellular settings
Blockage_real_cel = 32.4;
Blockage_est_cel = 32.4;
nu_real_cel = 2.1;
nu_est_cel = 2.1;
%noise parameters
if  nu_real_cel == 2.1
    sigma_cel = 4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_cel == 3.19
    sigma_cel = 7.82; % 4dB for LOS, 7.82 NLOS
end
%% RIS settings
Blockage_real_ris = 32.4;
Blockage_est_ris = 32.4;
nu_real_ris = 2.1;
nu_est_ris = 2.1;
%noise parameters
if  nu_real_ris == 2.1
    sigma_ris=4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_ris == 3.19
    sigma_ris=7.82; % 4dB for LOS, 7.82 NLOS
end
%% D2D settings
Blockage_real_d2d = 32.4;
Blockage_est_d2d = 32.4;
nu_real_d2d = 2.1;
nu_est_d2d = 2.1;
%noise parameters
if  nu_real_d2d == 2.1
    sigma_d2d = 4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_d2d == 3.19
    sigma_d2d = 7.82; % 4dB for LOS, 7.82 NLOS
end

%% CRLB
syms x
dLOS10 = double(solve(0.1 == 18/x + exp(-x/36)*(1 - 18/x), x));
dLOS20 = double(solve(0.2 == 18/x + exp(-x/36)*(1 - 18/x), x));
dLOS30 = double(solve(0.3 == 18/x + exp(-x/36)*(1 - 18/x), x));
dLOS70 = double(solve(0.7 == 18/x + exp(-x/36)*(1 - 18/x), x));

%% USER real position
index_counter = 1;
for realY = 0:10:100
%% add the cycle
    for realX = 0:10:100

        real_xy(index_counter,:) = [realX realY];
        %% 
        for index_mode = 1:length(mode_array)
            mode = mode_array(index_mode);
            if mode == 1 %RSS cellular

                X_a1 = 0;
                Y_a1 = 0;
                X_a2 = 0;
                Y_a2 = 100;
                X_a3 = 100;
                Y_a3 = 0;
                X_a4 = 100;
                Y_a4 = 100;
                X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];

                % Cellular
                %% CRLB from Fontanelli, Daniele, Farhad Shamsfakhr, and Luigi Palopoli. "Cramer–rao lower bound attainment in range-only positioning using geometry: The g-wls." IEEE Transactions on Instrumentation and Measurement 70 (2021): 1-14.
                % equation 13
                
                % pi_i is the distance between target node and anchor
                % x is the target node x coordonate
                % y is the target node y coordonate
                % X_i is the anchor x coordonate
                % Y_i is the anchor y coordonate
                
                % lambda_x1 = (x - X_i/pi_i)
                % lambda_y1 = (y - Y_i/pi_i)
                
                %% cellular CRLB
                pi_i1 = real(sqrt((X_a1- realX)^2 + (Y_a1- realY)^2)); %real distance between target node and anchor
                lambda_x1 = (realX - X_a1/pi_i1);
                lambda_y1 = (realY - Y_a1/pi_i1);
                
                pi_i2 = real(sqrt((X_a2- realX)^2 + (Y_a2- realY)^2)); %real distance between target node and anchor
                lambda_x2 = (realX - X_a2/pi_i2);
                lambda_y2 = (realY - Y_a2/pi_i2);
                
                pi_i3 = real(sqrt((X_a3- realX)^2 + (Y_a3- realY)^2)); %real distance between target node and anchor
                lambda_x3 = (realX - X_a3/pi_i3);
                lambda_y3 = (realY - Y_a3/pi_i3);

                pi_i4 = real(sqrt((X_a4- realX)^2 + (Y_a4- realY)^2)); %real distance between target node and anchor
                lambda_x4 = (realX - X_a4/pi_i4);
                lambda_y4 = (realY - Y_a4/pi_i4);
                
                H = [lambda_x1 lambda_y1; 
                    lambda_x2 lambda_y2;
                    lambda_x3 lambda_y3;
                    lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                
                % equation 12              
%                 nu_vector = [normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel)]; %zero-mean and Gaussian uncertainties
%                 nu = transpose(nu_vector);
%                 covariance_matrix = cov(nu.*nu',1); %% here is a mistake
                v = [sigma_cel sigma_cel  sigma_cel sigma_cel];
                covariance_matrix = diag(v);
                CRLB_cellular = inv(transpose(H)*(inv(covariance_matrix))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value(index_counter,:) = sqrt(CRLB_cellular(1,1) + CRLB_cellular(2,2));
                CRLB_cellular_xy(index_counter,:) = [CRLB_cellular(1,1) CRLB_cellular(2,2)];

                % CRLB with distance dependent noise (theoretical)
                d_array = [pi_i1 pi_i2 pi_i3 pi_i4];
                for i = 1:length(d_array)
                    d = d_array(i);
                    fc_GHz = 28;
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
                
                    % CRLB with distance dependent noise (experimental)
                    % empirical 10%
                    if d>dLOS10
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS10);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS10);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp10(i) = PL_LOS + PL_NLOS;
                    noise_var_emp10(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp10(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                
                    % empirical 20%
                    if d>dLOS20
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS20);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS20);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp20(i) = PL_LOS + PL_NLOS;
                    noise_var_emp20(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp20(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 30%
                    if d>dLOS30
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS30);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS30);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp30(i) = PL_LOS + PL_NLOS;
                    noise_var_emp30(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp30(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 70%
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
                
                covariance_matrix_dt = diag(noise_var);
                CRLB_cellular_dt = inv(transpose(H)*(inv(covariance_matrix_dt))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value_dt(index_counter,:) = sqrt(CRLB_cellular_dt(1,1) + CRLB_cellular_dt(2,2));
                CRLB_cellular_xy_dt(index_counter,:) = [CRLB_cellular_dt(1,1) CRLB_cellular_dt(2,2)];

                covariance_matrix_de10 = diag(noise_var_emp10);
                CRLB_cellular_de10 = inv(transpose(H)*(inv(covariance_matrix_de10))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value_de10(index_counter,:) = sqrt(CRLB_cellular_de10(1,1) + CRLB_cellular_de10(2,2));
                CRLB_cellular_xy_de10(index_counter,:) = [CRLB_cellular_de10(1,1) CRLB_cellular_de10(2,2)];

                covariance_matrix_de20 = diag(noise_var_emp20);
                CRLB_cellular_de20 = inv(transpose(H)*(inv(covariance_matrix_de20))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value_de20(index_counter,:) = sqrt(CRLB_cellular_de20(1,1) + CRLB_cellular_de20(2,2));
                CRLB_cellular_xy_de20(index_counter,:) = [CRLB_cellular_de20(1,1) CRLB_cellular_de20(2,2)];
                
                covariance_matrix_de30 = diag(noise_var_emp30);
                CRLB_cellular_de30 = inv(transpose(H)*(inv(covariance_matrix_de30))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value_de30(index_counter,:) = sqrt(CRLB_cellular_de30(1,1) + CRLB_cellular_de30(2,2));
                CRLB_cellular_xy_de30(index_counter,:) = [CRLB_cellular_de30(1,1) CRLB_cellular_de30(2,2)];

                covariance_matrix_de70 = diag(noise_var_emp70);
                CRLB_cellular_de70 = inv(transpose(H)*(inv(covariance_matrix_de70))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_cellular_value_de70(index_counter,:) = sqrt(CRLB_cellular_de70(1,1) + CRLB_cellular_de70(2,2));
                CRLB_cellular_xy_de70(index_counter,:) = [CRLB_cellular_de70(1,1) CRLB_cellular_de70(2,2)];
                
%                 GDoP_cellular = (inv(H'*H)); % equation 17
%                 GDoP_cellular_xy(index_counter,:) = [GDoP_cellular(1,1) GDoP_cellular(2,2)];
%                 GDoP_cellular_value(index_counter,:) = sqrt(GDoP_cellular(1,1) + GDoP_cellular(2,2));
% 
% %                    % equation 30 geometric lower bound:
%                 A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
%                      -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
%                      -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6
% 
%                 b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                      pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                      pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];
%  
%                 nu1 = normrnd(0,sigma_cel);
%                 nu2 = normrnd(0,sigma_cel);
%                 nu3 = normrnd(0,sigma_cel);
%                 nu4 = normrnd(0,sigma_cel);
%                 epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
%                            nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
%                            nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];
% 
%                 b_bar = b + epsilon; % equation 5
%                 eq_8 = [1+(sigma_cel*pi_i2^2)/(sigma_cel*pi_i1^2) 1 1
%                         1 1+(sigma_cel*pi_i3^2)/(sigma_cel*pi_i1^2) 1
%                         1 1 1+(sigma_cel*pi_i4^2)/(sigma_cel*pi_i1^2)];
%                 sum_epsilon = 4*sigma_cel*pi_i1^2*eq_8;
% %                    sum_epsilon = eq_8;
%                 p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
%                    
%                 pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x1 = (p_roof(1) - X_a1/pi_i1);
%                 lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
%                     
%                 pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x2 = (p_roof(1) - X_a2/pi_i2);
%                 lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
%                     
%                 pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x3 = (p_roof(1) - X_a3/pi_i3);
%                 lambda_y3 = (p_roof(2) - Y_a3/pi_i3);
%     
%                 pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x4 = (p_roof(1) - X_a4/pi_i4);
%                 lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
%                     
%                 H_bar = [lambda_x1 lambda_y1; 
%                 lambda_x2 lambda_y2;
%                 lambda_x3 lambda_y3;
%                 lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
%                     
%                 GLB_cel = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% % %                     GLB = sigma_cel*inv(H_bar*H_bar'); %equation 29
%                 GLB_cel_xy(index_counter,:) = [GLB_cel(1,1) GLB_cel(2,2)];
%                 GLB_cel_value(index_counter,:) = sqrt(GLB_cel(1,1) + GLB_cel(2,2));                
                
            elseif mode == 2 %RSS devices
                % user distribution in a circle of radius Rd
                if mode_d2d == 1
                %% random in all the area

                    X_B = X_BF;
                    Y_B = Y_BF;    
                    X = [X_B' Y_B']; 
                    X_a1 = X_B(1);
                    X_a2 = X_B(2);
                    X_a3 = X_B(3);
                    X_a4 = X_B(4);
                    Y_a1 = Y_B(1);
                    Y_a2 = Y_B(2);
                    Y_a3 = Y_B(3);
                    Y_a4 = Y_B(4);
                    
                elseif mode_d2d == 2

                    %% Random AROUND USER    
                    theta2 = rand(number_of_anchors,1)*2*pi;
                    X_B = realX + d0.*cos(theta2);
                    Y_B = realY + d0.*sin(theta2);

                    X_a1 = X_B(1);
                    X_a2 = X_B(2);
                    X_a3 = X_B(3);
                    X_a4 = X_B(4);
                    Y_a1 = Y_B(1);
                    Y_a2 = Y_B(2);
                    Y_a3 = Y_B(3);
                    Y_a4 = Y_B(4);
                    
                    X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];

                elseif mode_d2d==3
                  %%  fixed
                     
                    X_a1 = 7.1;
                    Y_a1 = 7.1;
                    X_a2 = 7.1;
                    Y_a2 = 92.9;
                    X_a3 = 92.9;
                    Y_a3 = 7.1;
                    X_a4 = 92.9;
                    Y_a4 = 92.9;
                    
                    X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];
 
                elseif mode_d2d==4
                   %% corner variable 
                    X_a1 = realX - randi([5 20]);
                    Y_a1 = realY - randi([5 20]);
                    X_a2 = realX + randi([5 20]);
                    Y_a2 = realY - randi([5 20]);
                    X_a3 = realX - randi([5 20]);
                    Y_a3 = realY + randi([5 20]);
                    X_a4 = realX + randi([5 20]);
                    Y_a4 = realY + randi([5 20]);

                    X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];
                end

                % D2D
                %% CRLB from Fontanelli, Daniele, Farhad Shamsfakhr, and Luigi Palopoli. "Cramer–rao lower bound attainment in range-only positioning using geometry: The g-wls." IEEE Transactions on Instrumentation and Measurement 70 (2021): 1-14.
                % equation 13
                    
                % pi_i is the distance between target node and anchor
                % x is the target node x coordonate
                % y is the target node y coordonate
                % X_i is the anchor x coordonate
                % Y_i is the anchor y coordonate
                    
                % lambda_x1 = (x - X_i/pi_i)
                % lambda_y1 = (y - Y_i/pi_i)
                    
                %% D2D CRLB
                pi_i1 = real(sqrt((X_a1- realX)^2 + (Y_a1- realY)^2)); %real distance between target node and anchor
                lambda_x1 = (realX - X_a1/pi_i1);
                lambda_y1 = (realY - Y_a1/pi_i1);
                    
                pi_i2 = real(sqrt((X_a2- realX)^2 + (Y_a2- realY)^2)); %real distance between target node and anchor
                lambda_x2 = (realX - X_a2/pi_i2);
                lambda_y2 = (realY - Y_a2/pi_i2);
                    
                pi_i3 = real(sqrt((X_a3- realX)^2 + (Y_a3- realY)^2)); %real distance between target node and anchor
                lambda_x3 = (realX - X_a3/pi_i3);
                lambda_y3 = (realY - Y_a3/pi_i3);
    
                pi_i4 = real(sqrt((X_a4- realX)^2 + (Y_a4- realY)^2)); %real distance between target node and anchor
                lambda_x4 = (realX - X_a4/pi_i4);
                lambda_y4 = (realY - Y_a4/pi_i4);
                    
                H = [lambda_x1 lambda_y1; 
                     lambda_x2 lambda_y2;
                     lambda_x3 lambda_y3;
                     lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                    
                % equation 12              
    %                 nu_vector = [normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel)]; %zero-mean and Gaussian uncertainties
    %                 nu = transpose(nu_vector);
    %                 covariance_matrix = cov(nu.*nu',1); 
                v = [sigma_d2d sigma_d2d  sigma_d2d sigma_d2d];
                covariance_matrix = diag(v);
                CRLB_d2d = inv(transpose(H)*(inv(covariance_matrix))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value(index_counter,:) = sqrt(CRLB_d2d(1,1) + CRLB_d2d(2,2));
                CRLB_d2d_xy(index_counter,:) = [CRLB_d2d(1,1) CRLB_d2d(2,2)];

                % CRLB with distance dependent noise (theoretical)
                d_array = [pi_i1 pi_i2 pi_i3 pi_i4];
                for i = 1:length(d_array)
                    d = d_array(i);
                    fc_GHz = 28;
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
                    Pt = db2pow(10)/1000;
                    noise_var(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_theor(i)) / ( (100*10^6)^2 *(Pt/Pn) ) );
                
                    % CRLB with distance dependent noise (experimental)
                    % empirical 10%
                    if d>dLOS10
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS10);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS10);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp10(i) = PL_LOS + PL_NLOS;
                    noise_var_emp10(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp10(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                
                    % empirical 20%
                    if d>dLOS20
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS20);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS20);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp20(i) = PL_LOS + PL_NLOS;
                    noise_var_emp20(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp20(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 30%
                    if d>dLOS30
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS30);
                        PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS30);
                    else
                        PL_NLOS = 0;
                        PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);
                    end
                    PL_emp30(i) = PL_LOS + PL_NLOS;
                    noise_var_emp30(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp30(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 70%
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

                covariance_matrix_dt = diag(noise_var);
                CRLB_d2d_dt = inv(transpose(H)*(inv(covariance_matrix_dt))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value_dt(index_counter,:) = sqrt(CRLB_d2d_dt(1,1) + CRLB_d2d_dt(2,2));
                CRLB_d2d_xy_dt(index_counter,:) = [CRLB_d2d_dt(1,1) CRLB_d2d_dt(2,2)];

                covariance_matrix_de10 = diag(noise_var_emp10);
                CRLB_d2d_de10 = inv(transpose(H)*(inv(covariance_matrix_de10))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value_de10(index_counter,:) = sqrt(CRLB_d2d_de10(1,1) + CRLB_d2d_de10(2,2));
                CRLB_d2d_xy_de10(index_counter,:) = [CRLB_d2d_de10(1,1) CRLB_d2d_de10(2,2)];
                
                covariance_matrix_de20 = diag(noise_var_emp20);
                CRLB_d2d_de20 = inv(transpose(H)*(inv(covariance_matrix_de20))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value_de20(index_counter,:) = sqrt(CRLB_d2d_de20(1,1) + CRLB_d2d_de20(2,2));
                CRLB_d2d_xy_de20(index_counter,:) = [CRLB_d2d_de20(1,1) CRLB_d2d_de20(2,2)];

                covariance_matrix_de30 = diag(noise_var_emp30);
                CRLB_d2d_de30 = inv(transpose(H)*(inv(covariance_matrix_de30))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value_de30(index_counter,:) = sqrt(CRLB_d2d_de30(1,1) + CRLB_d2d_de30(2,2));
                CRLB_d2d_xy_de30(index_counter,:) = [CRLB_d2d_de30(1,1) CRLB_d2d_de30(2,2)];

                covariance_matrix_de70 = diag(noise_var_emp70);
                CRLB_d2d_de70 = inv(transpose(H)*(inv(covariance_matrix_de70))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_d2d_value_de70(index_counter,:) = sqrt(CRLB_d2d_de70(1,1) + CRLB_d2d_de70(2,2));
                CRLB_d2d_xy_de70(index_counter,:) = [CRLB_d2d_de70(1,1) CRLB_d2d_de70(2,2)];
                
%                 GDoP_d2d = (inv(H'*H)); % equation 17
%                 GDoP_d2d_xy(index_counter,:) = [GDoP_d2d(1,1) GDoP_d2d(2,2)];
%                 GDoP_d2d_value(index_counter,:) = sqrt(GDoP_d2d(1,1) + GDoP_d2d(2,2));
% % 
% %               % equation 30 geometric lower bound:
%                 A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
%                     -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
%                     -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6
% 
%                 b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                     pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                     pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];
%  
%                 nu1 = normrnd(0,sigma_d2d);
%                 nu2 = normrnd(0,sigma_d2d);
%                 nu3 = normrnd(0,sigma_d2d);
%                 nu4 = normrnd(0,sigma_d2d);
%                 epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
%                           nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
%                           nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];
% 
%                 b_bar = b + epsilon; % equation 5
%                 eq_8 = [1+(sigma_d2d*pi_i2^2)/(sigma_d2d*pi_i1^2) 1 1
%                        1 1+(sigma_d2d*pi_i3^2)/(sigma_d2d*pi_i1^2) 1
%                        1 1 1+(sigma_d2d*pi_i4^2)/(sigma_d2d*pi_i1^2)];
%                 sum_epsilon = 4*sigma_d2d*pi_i1^2*eq_8;
% %                    sum_epsilon = eq_8;
%                 p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
%                    
%                 pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x1 = (p_roof(1) - X_a1/pi_i1);
%                 lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
%                 
%                 pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x2 = (p_roof(1) - X_a2/pi_i2);
%                 lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
%                 
%                 pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x3 = (p_roof(1) - X_a3/pi_i3);
%                 lambda_y3 = (p_roof(2) - Y_a3/pi_i3);
% 
%                 pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x4 = (p_roof(1) - X_a4/pi_i4);
%                 lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
%                     
%                 H_bar = [lambda_x1 lambda_y1; 
%                     lambda_x2 lambda_y2;
%                     lambda_x3 lambda_y3;
%                     lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
%                 
%                 GLB_d2d = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% % %                     GLB = sigma_cel*inv(H_bar*H_bar'); %equation 29
%                 GLB_d2d_xy(index_counter,:) = [GLB_d2d(1,1) GLB_d2d(2,2)];
%                 GLB_d2d_value(index_counter,:) = sqrt(GLB_d2d(1,1) + GLB_d2d(2,2));

             elseif mode == 3 %% RISs
                % RISs are located clode to the BSs
                X_a1 = 7.1;
                Y_a1 = 7.1;
                X_a2 = 7.1;
                Y_a2 = 92.9;
                X_a3 = 92.9;
                Y_a3 = 7.1;
                X_a4 = 92.9;
                Y_a4 = 92.9;

                X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];
            %     plot (X_a1,Y_a1,'o',X_a2,Y_a2,'*',X_a3,Y_a3,'s',X_a4,Y_a4,'^')

                %BS positions
                X_a1_BS = 0;
                Y_a1_BS = 0;
                X_a2_BS = 0;
                Y_a2_BS = 100;
                X_a3_BS = 100;
                Y_a3_BS = 0;
                X_a4_BS = 100;
                Y_a4_BS = 100;

                % RIS
                %% CRLB from Fontanelli, Daniele, Farhad Shamsfakhr, and Luigi Palopoli. "Cramer–rao lower bound attainment in range-only positioning using geometry: The g-wls." IEEE Transactions on Instrumentation and Measurement 70 (2021): 1-14.
                % equation 13
                
                % pi_i is the distance between target node and anchor
                % x is the target node x coordonate
                % y is the target node y coordonate
                % X_i is the anchor x coordonate
                % Y_i is the anchor y coordonate
                
                % lambda_x1 = (x - X_i/pi_i)
                % lambda_y1 = (y - Y_i/pi_i)
                
                %% RIS CRLB
                pi_i1 = real(sqrt((X_a1- realX)^2 + (Y_a1- realY)^2)); %real distance between target node and anchor
                lambda_x1 = (realX - X_a1/pi_i1);
                lambda_y1 = (realY - Y_a1/pi_i1);
                
                pi_i2 = real(sqrt((X_a2- realX)^2 + (Y_a2- realY)^2)); %real distance between target node and anchor
                lambda_x2 = (realX - X_a2/pi_i2);
                lambda_y2 = (realY - Y_a2/pi_i2);
                
                pi_i3 = real(sqrt((X_a3- realX)^2 + (Y_a3- realY)^2)); %real distance between target node and anchor
                lambda_x3 = (realX - X_a3/pi_i3);
                lambda_y3 = (realY - Y_a3/pi_i3);

                pi_i4 = real(sqrt((X_a4- realX)^2 + (Y_a4- realY)^2)); %real distance between target node and anchor
                lambda_x4 = (realX - X_a4/pi_i4);
                lambda_y4 = (realY - Y_a4/pi_i4);
                
                H = [lambda_x1 lambda_y1; 
                    lambda_x2 lambda_y2;
                    lambda_x3 lambda_y3;
                    lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                
                % equation 12              
%                 nu_vector = [normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel) normrnd(0,sigma_cel)]; %zero-mean and Gaussian uncertainties
%                 nu = transpose(nu_vector);
%                 covariance_matrix = cov(nu.*nu',1); 
                v = [sigma_ris sigma_ris  sigma_ris sigma_ris];
                covariance_matrix = diag(v);
                CRLB_RIS = inv(transpose(H)*(inv(covariance_matrix))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value(index_counter,:) = sqrt(CRLB_RIS(1,1) + CRLB_RIS(2,2));
                CRLB_RIS_xy(index_counter,:) = [CRLB_RIS(1,1) CRLB_RIS(2,2)];

                % CRLB with distance dependent noise (theoretical)
                d_array = [pi_i1 pi_i2 pi_i3 pi_i4];
                for i = 1:length(d_array)
                    d = d_array(i);
                    fc_GHz = 28;
%                     FSPL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d); % dB 
%                     FSPL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d); % dB 

                    FSPL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d^(2.1)); % linear n
                    FSPL_NLOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d)^(3.19); % linear 
                    FSPL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                    FSPL_LOS = pow2db(((sqrt((1/(FSPL_LOS_SR*FSPL_LOS_RD))))*1024)^(-2));
                    FSPL_NLOS = pow2db(((sqrt((1/(FSPL_LOS_SR*FSPL_NLOS_RD))))*1024)^(-2));
                    
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
                
                    % CRLB with distance dependent noise (experimental)
                    % empirical 10%
                    if d>dLOS10
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS10);
%                         PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS10);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(dLOS10^(2.1)); % linear n
                        PL_NLOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d-dLOS10)^(3.19); % linear 
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_NLOS_RD))))*1024)^(-2));

                    else
%                         PL_NLOS = 0;
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d^(2.1)); % linear n
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = 0;


                    end
                    PL_emp10(i) = PL_LOS + PL_NLOS;
                    noise_var_emp10(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp10(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                
                    % empirical 20%
                    if d>dLOS20
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS20);
%                         PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS20);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(dLOS20^(2.1)); % linear n
                        PL_NLOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d-dLOS20)^(3.19); % linear 
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_NLOS_RD))))*1024)^(-2));

                    else
%                         PL_NLOS = 0;
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d^(2.1)); % linear n
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = 0;

                    end
                    PL_emp20(i) = PL_LOS + PL_NLOS;
                    noise_var_emp20(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp20(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 30%
                    if d>dLOS30
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS30);
%                         PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS30);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(dLOS30^(2.1)); % linear n
                        PL_NLOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d-dLOS30)^(3.19); % linear 
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_NLOS_RD))))*1024)^(-2));
                    else
%                         PL_NLOS = 0;
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d^(2.1)); % linear n
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = 0;
                    end
                    PL_emp30(i) = PL_LOS + PL_NLOS;
                    noise_var_emp30(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp30(i)) / ( (100*10^6)^2 *(Pt/Pn) ));
                    % empirical 70%
                    if d>dLOS70
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(dLOS70);
%                         PL_NLOS = 32.4 + 20*log10(fc_GHz) + 31.9*log10(d-dLOS70);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(dLOS70^(2.1)); % linear n
                        PL_NLOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d-dLOS70)^(3.19); % linear 
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_NLOS_RD))))*1024)^(-2));
                    else
%                         PL_NLOS = 0;
%                         PL_LOS = 32.4 + 20*log10(fc_GHz) + 21*log10(d);

                        PL_LOS_RD = 10^(2*log10(fc_GHz) + 3.24)*(d^(2.1)); % linear n
                        PL_LOS_SR = 10^(2*log10(fc_GHz) + 3.24)*(10^(2.1)); % linear n
                    
                        PL_LOS = pow2db(((sqrt((1/(PL_LOS_SR*PL_LOS_RD))))*1024)^(-2));
                        PL_NLOS = 0;
                    end
                    PL_emp70(i) = PL_LOS + PL_NLOS;
                    noise_var_emp70(i) = pow2db(physconst('LightSpeed')^2*db2pow(PL_emp70(i)) / ( (100*10^6)^2 *(Pt/Pn) ));

                end
                
                covariance_matrix_dt = diag(noise_var);
                CRLB_RIS_dt = inv(transpose(H)*(inv(covariance_matrix_dt))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value_dt(index_counter,:) = sqrt(CRLB_RIS_dt(1,1) + CRLB_RIS_dt(2,2));
                CRLB_RIS_xy_dt(index_counter,:) = [CRLB_RIS_dt(1,1) CRLB_RIS_dt(2,2)];

                covariance_matrix_de10 = diag(noise_var_emp10);
                CRLB_RIS_de10 = inv(transpose(H)*(inv(covariance_matrix_de10))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value_de10(index_counter,:) = sqrt(CRLB_RIS_de10(1,1) + CRLB_RIS_de10(2,2));
                CRLB_RIS_xy_de10(index_counter,:) = [CRLB_RIS_de10(1,1) CRLB_RIS_de10(2,2)];

                covariance_matrix_de20 = diag(noise_var_emp20);
                CRLB_RIS_de20 = inv(transpose(H)*(inv(covariance_matrix_de20))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value_de20(index_counter,:) = sqrt(CRLB_RIS_de20(1,1) + CRLB_RIS_de20(2,2));
                CRLB_RIS_xy_de20(index_counter,:) = [CRLB_RIS_de20(1,1) CRLB_RIS_de20(2,2)];

                covariance_matrix_de30 = diag(noise_var_emp30);
                CRLB_RIS_de30 = inv(transpose(H)*(inv(covariance_matrix_de30))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value_de30(index_counter,:) = sqrt(CRLB_RIS_de30(1,1) + CRLB_RIS_de30(2,2));
                CRLB_RIS_xy_de30(index_counter,:) = [CRLB_RIS_de30(1,1) CRLB_RIS_de30(2,2)];

                covariance_matrix_de70 = diag(noise_var_emp70);
                CRLB_RIS_de70 = inv(transpose(H)*(inv(covariance_matrix_de70))*H);
                % output in the diagonal is the error on x and y axis
                CRLB_RIS_value_de70(index_counter,:) = sqrt(CRLB_RIS_de70(1,1) + CRLB_RIS_de70(2,2));
                CRLB_RIS_xy_de70(index_counter,:) = [CRLB_RIS_de70(1,1) CRLB_RIS_de70(2,2)];
                    
%                 GDoP_RIS = (inv(H'*H)); % equation 17
%                 GDoP_RIS_xy(index_counter,:) = [GDoP_RIS(1,1) GDoP_RIS(2,2)];
%                 GDoP_RIS_value(index_counter,:) = sqrt(GDoP_RIS(1,1) + GDoP_RIS(2,2));
% 
% %                    % equation 30 geometric lower bound:
%                 A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
%                      -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
%                      -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6
% 
%                 b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                      pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
%                      pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];
%  
%                 nu1 = normrnd(0,sigma_ris);
%                 nu2 = normrnd(0,sigma_ris);
%                 nu3 = normrnd(0,sigma_ris);
%                 nu4 = normrnd(0,sigma_ris);
%                 epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
%                            nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
%                            nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];
% 
%                 b_bar = b + epsilon; % equation 5
%                 eq_8 = [1+(sigma_ris*pi_i2^2)/(sigma_ris*pi_i1^2) 1 1
%                         1 1+(sigma_ris*pi_i3^2)/(sigma_ris*pi_i1^2) 1
%                         1 1 1+(sigma_ris*pi_i4^2)/(sigma_ris*pi_i1^2)];
%                 sum_epsilon = 4*sigma_ris*pi_i1^2*eq_8;
% %                    sum_epsilon = eq_8;
%                 p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
%                    
%                 pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x1 = (p_roof(1) - X_a1/pi_i1);
%                 lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
%                     
%                 pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x2 = (p_roof(1) - X_a2/pi_i2);
%                 lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
%                     
%                 pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x3 = (p_roof(1) - X_a3/pi_i3);
%                 lambda_y3 = (p_roof(2) - Y_a3/pi_i3);
%     
%                 pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
%                 lambda_x4 = (p_roof(1) - X_a4/pi_i4);
%                 lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
%                     
%                 H_bar = [lambda_x1 lambda_y1; 
%                         lambda_x2 lambda_y2;
%                         lambda_x3 lambda_y3;
%                         lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
%                     
%                 GLB_RIS = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% % %                     GLB = sigma_ris*inv(H_bar*H_bar'); %equation 29
%                 GLB_RIS_xy(index_counter,:) = [GLB_RIS(1,1) GLB_RIS(2,2)];
%                 GLB_RIS_value(index_counter,:) = sqrt(GLB_RIS(1,1) + GLB_RIS(2,2));                 
            
            end
        end
        index_counter=index_counter+1;
    end
      
end      
%real
CRLB_cellular_value_de10 = real(CRLB_cellular_value_de10);
CRLB_cellular_value_de20 = real(CRLB_cellular_value_de20);
CRLB_cellular_value_de30 = real(CRLB_cellular_value_de30);
CRLB_cellular_value_de70 = real(CRLB_cellular_value_de70);
CRLB_cellular_value_dt = real(CRLB_cellular_value_dt);

CRLB_d2d_value_de10 = real(CRLB_d2d_value_de10);
CRLB_d2d_value_de20 = real(CRLB_d2d_value_de20);
CRLB_d2d_value_de30 = real(CRLB_d2d_value_de30);
CRLB_d2d_value_de70 = real(CRLB_d2d_value_de70);
CRLB_d2d_value_dt = real(CRLB_d2d_value_dt);

CRLB_RIS_value_de10 = real(CRLB_RIS_value_de10);
CRLB_RIS_value_de20 = real(CRLB_RIS_value_de20);
CRLB_RIS_value_de30 = real(CRLB_RIS_value_de30);
CRLB_RIS_value_de70 = real(CRLB_RIS_value_de70);
CRLB_RIS_value_dt = real(CRLB_RIS_value_dt);


% 50th 75th 95th percentile
% lower bounds
error_CRLB_cel_prc = prctile(CRLB_cellular_value ,[50 75 95]);
error_CRLB_d2d_prc = prctile(CRLB_d2d_value ,[50 75 95]);
error_CRLB_RIS_prc = prctile(CRLB_RIS_value ,[50 75 95]);
error_CRLB_cel_prc_dt = prctile(CRLB_cellular_value_dt ,[50 75 95]);
error_CRLB_cel_prc_de10 = prctile(CRLB_cellular_value_de10 ,[50 75 95]);
error_CRLB_cel_prc_de20 = prctile(CRLB_cellular_value_de20 ,[50 75 95]);
error_CRLB_cel_prc_de30 = prctile(CRLB_cellular_value_de30 ,[50 75 95]);
error_CRLB_cel_prc_de70 = prctile(CRLB_cellular_value_de70 ,[50 75 95]);

error_CRLB_d2d_prc_dt = prctile(CRLB_d2d_value_dt ,[50 75 95]);
error_CRLB_d2d_prc_de10 = prctile(CRLB_d2d_value_de10 ,[50 75 95]);
error_CRLB_d2d_prc_de20 = prctile(CRLB_d2d_value_de20 ,[50 75 95]);
error_CRLB_d2d_prc_de30 = prctile(CRLB_d2d_value_de30 ,[50 75 95]);
error_CRLB_d2d_prc_de70 = prctile(CRLB_d2d_value_de70 ,[50 75 95]);

error_CRLB_RIS_prc_dt = prctile(CRLB_RIS_value_dt ,[50 75 95]);
error_CRLB_RIS_prc_de10 = prctile(CRLB_RIS_value_de10 ,[50 75 95]);
error_CRLB_RIS_prc_de20 = prctile(CRLB_RIS_value_de20 ,[50 75 95]);
error_CRLB_RIS_prc_de30 = prctile(CRLB_RIS_value_de30 ,[50 75 95]);
error_CRLB_RIS_prc_de70 = prctile(CRLB_RIS_value_de70 ,[50 75 95]);
% error_GDoP_cel_prc = prctile(GDoP_cellular_value ,[50 75 95]);
% error_GDoP_d2d_prc = prctile(GDoP_d2d_value ,[50 75 95]);
% error_GDoP_RIS_prc = prctile(GDoP_RIS_value ,[50 75 95]);
% error_GLB_cel_prc = prctile(GLB_cel_value ,[50 75 95]);
% error_GLB_d2d_prc = prctile(GLB_d2d_value ,[50 75 95]);
% error_GLB_RIS_prc = prctile(GLB_RIS_value ,[50 75 95]);


%MSE 
error_CRLB_cel_MSE = mean(CRLB_cellular_value.^2);
error_CRLB_d2d_MSE = mean(CRLB_d2d_value.^2);
error_CRLB_RIS_MSE = mean(CRLB_RIS_value.^2);
error_CRLB_cel_MSE_dt = mean(CRLB_cellular_value_dt.^2);
error_CRLB_cel_MSE_de10 = mean(CRLB_cellular_value_de10.^2);
error_CRLB_cel_MSE_de20 = mean(CRLB_cellular_value_de20.^2);
error_CRLB_cel_MSE_de30 = mean(CRLB_cellular_value_de30.^2);
error_CRLB_cel_MSE_de70 = mean(CRLB_cellular_value_de70.^2);
error_CRLB_d2d_MSE_dt = mean(CRLB_d2d_value_dt.^2);
error_CRLB_d2d_MSE_de10 = mean(CRLB_d2d_value_de10.^2);
error_CRLB_d2d_MSE_de20 = mean(CRLB_d2d_value_de20.^2);
error_CRLB_d2d_MSE_de30 = mean(CRLB_d2d_value_de30.^2);
error_CRLB_d2d_MSE_de70 = mean(CRLB_d2d_value_de70.^2);

error_CRLB_RIS_MSE_dt = mean(CRLB_RIS_value_dt.^2);
error_CRLB_RIS_MSE_de10 = mean(CRLB_RIS_value_de10.^2);
error_CRLB_RIS_MSE_de20 = mean(CRLB_RIS_value_de20.^2);
error_CRLB_RIS_MSE_de30 = mean(CRLB_RIS_value_de30.^2);
error_CRLB_RIS_MSE_de70 = mean(CRLB_RIS_value_de70.^2);
% error_GDoP_cel_MSE = mean(GDoP_cellular_value.^2);
% error_GDoP_d2d_MSE = mean(GDoP_d2d_value.^2);
% error_GDoP_RIS_MSE = mean(GDoP_RIS_value.^2);
% error_GLB_cel_MSE = mean(GLB_cel_value.^2);
% error_GLB_d2d_MSE = mean(GLB_d2d_value.^2);
% error_GLB_RIS_MSE = mean(GLB_RIS_value.^2);

%RMSE
error_CRLB_cel_RMSE = sqrt(error_CRLB_cel_MSE);
error_CRLB_d2d_RMSE = sqrt(error_CRLB_d2d_MSE);
error_CRLB_RIS_RMSE = sqrt(error_CRLB_RIS_MSE);
error_CRLB_cel_RMSE_dt = sqrt(error_CRLB_cel_MSE_dt);
error_CRLB_cel_RMSE_de10 = sqrt(error_CRLB_cel_MSE_de10);
error_CRLB_cel_RMSE_de20 = sqrt(error_CRLB_cel_MSE_de20);
error_CRLB_cel_RMSE_de30 = sqrt(error_CRLB_cel_MSE_de30);
error_CRLB_cel_RMSE_de70 = sqrt(error_CRLB_cel_MSE_de70);

error_CRLB_d2d_RMSE_dt = sqrt(error_CRLB_d2d_MSE_dt);
error_CRLB_d2d_RMSE_de10 = sqrt(error_CRLB_d2d_MSE_de10);
error_CRLB_d2d_RMSE_de20 = sqrt(error_CRLB_d2d_MSE_de20);
error_CRLB_d2d_RMSE_de30 = sqrt(error_CRLB_d2d_MSE_de30);
error_CRLB_d2d_RMSE_de70 = sqrt(error_CRLB_d2d_MSE_de70);

error_CRLB_RIS_RMSE_dt = sqrt(error_CRLB_RIS_MSE_dt);
error_CRLB_RIS_RMSE_de10 = sqrt(error_CRLB_RIS_MSE_de10);
error_CRLB_RIS_RMSE_de20 = sqrt(error_CRLB_RIS_MSE_de20);
error_CRLB_RIS_RMSE_de30 = sqrt(error_CRLB_RIS_MSE_de30);
error_CRLB_RIS_RMSE_de70 = sqrt(error_CRLB_RIS_MSE_de70);

% error_GDoP_cel_RMSE = sqrt(error_GDoP_cel_MSE);
% error_GDoP_d2d_RMSE = sqrt(error_GDoP_d2d_MSE);
% error_GDoP_RIS_RMSE = sqrt(error_GDoP_RIS_MSE);
% error_GLB_cel_RMSE = sqrt(error_GLB_cel_MSE);
% error_GLB_d2d_RMSE = sqrt(error_GLB_d2d_MSE);
% error_GLB_RIS_RMSE = sqrt(error_GLB_RIS_MSE);

set(0,'DefaultAxesFontSize', 13)
% Change default text fonts.
set(0,'DefaultTextFontSize', 13)
%CRLB
figure(1)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_reshaped = reshape(CRLB_cellular_value, [11, 11]);
contour(X,Y,CRLB_cellular_value_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(2)
CRLB_d2d_value_reshaped = reshape(CRLB_d2d_value, [11, 11]);
contour(X,Y,CRLB_d2d_value_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(3)
CRLB_RIS_value_reshaped = reshape(CRLB_RIS_value, [11, 11]);
contour(X,Y,CRLB_RIS_value_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% CRLB distance-dependent theoretical
figure(4)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_dt_reshaped = reshape(CRLB_cellular_value_dt, [11, 11]);
contour(X,Y,CRLB_cellular_value_dt_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(5)
CRLB_d2d_value_dt_reshaped = reshape(CRLB_d2d_value_dt, [11, 11]);
contour(X,Y,CRLB_d2d_value_dt_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(6)
CRLB_RIS_value_dt_reshaped = reshape(CRLB_RIS_value_dt, [11, 11]);
contour(X,Y,CRLB_RIS_value_dt_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% CRLB distance-dependent 10%
figure(7)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_de10_reshaped = reshape(CRLB_cellular_value_de10, [11, 11]);
contour(X,Y,CRLB_cellular_value_de10_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(8)
CRLB_d2d_value_de10_reshaped = reshape(CRLB_d2d_value_de10, [11, 11]);
contour(X,Y,CRLB_d2d_value_de10_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(9)
CRLB_RIS_value_de10_reshaped = reshape(CRLB_RIS_value_de10, [11, 11]);
contour(X,Y,CRLB_RIS_value_de10_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% CRLB distance-dependent 20%
figure(10)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_de20_reshaped = reshape(CRLB_cellular_value_de20, [11, 11]);
contour(X,Y,CRLB_cellular_value_de20_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(11)
CRLB_d2d_value_de20_reshaped = reshape(CRLB_d2d_value_de20, [11, 11]);
contour(X,Y,CRLB_d2d_value_de20_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(12)
CRLB_RIS_value_de20_reshaped = reshape(CRLB_RIS_value_de20, [11, 11]);
contour(X,Y,CRLB_RIS_value_de20_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% CRLB distance-dependent 30%
figure(13)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_de30_reshaped = reshape(CRLB_cellular_value_de30, [11, 11]);
contour(X,Y,CRLB_cellular_value_de30_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(14)
CRLB_d2d_value_de30_reshaped = reshape(CRLB_d2d_value_de30, [11, 11]);
contour(X,Y,CRLB_d2d_value_de30_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(15)
CRLB_RIS_value_de30_reshaped = reshape(CRLB_RIS_value_de30, [11, 11]);
contour(X,Y,CRLB_RIS_value_de30_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% CRLB distance-dependent 70%
figure(16)
[X,Y] = meshgrid(linspace(0,100, 11),linspace(0,100, 11));
CRLB_cellular_value_de70_reshaped = reshape(CRLB_cellular_value_de70, [11, 11]);
contour(X,Y,CRLB_cellular_value_de70_reshaped,'LineWidth', 3, 'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(17)
CRLB_d2d_value_de70_reshaped = reshape(CRLB_d2d_value_de70, [11, 11]);
contour(X,Y,CRLB_d2d_value_de70_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")
figure(18)
CRLB_RIS_value_de70_reshaped = reshape(CRLB_RIS_value_de70, [11, 11]);
contour(X,Y,CRLB_RIS_value_de70_reshaped,'LineWidth', 3,'ShowText','on')
xlabel('X, m')
ylabel('Y, m')
xticks(linspace(0,100, 11))
yticks(linspace(0,100, 11))
axis("equal")

% % GDoP
% figure(4)
% GDoP_cellular_value_reshaped = reshape(GDoP_cellular_value, [11, 11]);
% contour(X,Y,GDoP_cellular_value_reshaped,'LineWidth', 3, 'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")
% figure(5)
% GDoP_d2d_value_reshaped = reshape(GDoP_d2d_value, [11, 11]);
% contour(X,Y,GDoP_d2d_value_reshaped,'LineWidth', 3,'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")
% figure(6)
% GDoP_RIS_value_reshaped = reshape(GDoP_RIS_value, [11, 11]);
% contour(X,Y,GDoP_RIS_value_reshaped,'LineWidth', 3,'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")
% % GLB
% figure(7)
% GLB_cellular_value_reshaped = reshape(GLB_cel_value, [11, 11]);
% contour(X,Y,GLB_cellular_value_reshaped,'LineWidth', 3, 'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")
% figure(8)
% GLB_d2d_value_reshaped = reshape(GLB_d2d_value, [11, 11]);
% contour(X,Y,GLB_d2d_value_reshaped,'LineWidth', 3,'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")
% figure(9)
% GLB_RIS_value_reshaped = reshape(GLB_RIS_value, [11, 311]);
% contour(X,Y,GLB_RIS_value_reshaped,'LineWidth', 3,'ShowText','on')
% xlabel('X, m')
% ylabel('Y, m')
% xticks(linspace(0,100, 11))
% yticks(linspace(0,100, 11))
% axis("equal")

% save('figure5_new.mat'); 
% load('figure5.mat')