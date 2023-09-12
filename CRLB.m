clear all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Table IV Lower bounds, [m]. Saved data are named as bounds_d2d....                                                        %
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 3.0 (Last edited: 2023-09-06)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RSS  location estimation
%% parameters
prompt_number_of_evaluation_points = "Enter number of evaluation points (M): ";
number_of_evaluation_points = input(prompt_number_of_evaluation_points);
%number_of_evaluation_points = 30;
prompt_runs = "Enter number of experiment repetitions (rep): ";
runs = input(prompt_runs);
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
NN=1024; % Number of reflective elements
%% D2D mode switching
prompt_d2d_mode = "Choose D2D user distribution? Options are 1,2,3,4: ";
mode_d2d = input(prompt_d2d_mode);
%mode_d2d = 1;
%1 - random 4 points within the considered area
%2 - random 4 points near real point within the radius of 20 m
d0 = 20; %radius of 20 m
%3 - equal placement of RISs and relays
%4 - 4 points on the corners on variable (5,20) distance [not used in the
%paper]
%% Cellular settings
prompt_Blockage_real_cel = "CELLULAR setup: Enter real non-blocked/blocked path (32.4/47.4): ";
Blockage_real_cel = input(prompt_Blockage_real_cel);
prompt_Blockage_est_cel = "CELLULAR setup: Enter estimated non-blocked/blocked path (32.4/47.4): ";
Blockage_est_cel = input(prompt_Blockage_est_cel);
prompt_nu_real_cel = "CELLULAR setup: Enter real PLE LOS/NLOS (2.1/3.19): ";
nu_real_cel = input(prompt_nu_real_cel);
prompt_nu_est_cel = "CELLULAR setup: Enter estimated PLE LOS/NLOS (2.1/3.19): ";
nu_est_cel = input(prompt_nu_est_cel);
%noise parameters
if  nu_real_cel == 2.1
    sigma_cel = 4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_cel == 3.19
    sigma_cel = 7.82; % 4dB for LOS, 7.82 NLOS
end
%% RIS settings
prompt_Blockage_real_ris = "RIS setup: Enter real non-blocked/blocked path (32.4/47.4): ";
Blockage_real_ris = input(prompt_Blockage_real_ris);
prompt_Blockage_est_ris = "RIS setup: Enter estimated non-blocked/blocked path (32.4/47.4): ";
Blockage_est_ris = input(prompt_Blockage_est_ris);
prompt_nu_real_ris = "RIS setup: Enter real PLE LOS/NLOS (2.1/3.19): ";
nu_real_ris = input(prompt_nu_real_ris);
prompt_nu_est_ris = "RIS setup: Enter estimated PLE LOS/NLOS (2.1/3.19): ";
nu_est_ris = input(prompt_nu_est_ris);
%noise parameters
if  nu_real_ris == 2.1
    sigma_ris = 4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_ris == 3.19
    sigma_ris = 7.82; % 4dB for LOS, 7.82 NLOS
end
%% D2D settings
prompt_Blockage_real_d2d = "D2D setup: Enter real non-blocked/blocked path (32.4/47.4): ";
Blockage_real_d2d = input(prompt_Blockage_real_d2d);
prompt_Blockage_est_d2d = "D2D setup: Enter estimated non-blocked/blocked path (32.4/47.4): ";
Blockage_est_d2d = input(prompt_Blockage_est_d2d);
prompt_nu_real_d2d = "D2D setup: Enter real PLE LOS/NLOS (2.1/3.19): ";
nu_real_d2d = input(prompt_nu_real_d2d);
prompt_nu_est_d2d = "D2D setup: Enter estimated PLE LOS/NLOS (2.1/3.19): ";
nu_est_d2d = input(prompt_nu_est_d2d);
%noise parameters
if  nu_real_d2d == 2.1
    sigma_d2d = 4; % 4dB for LOS, 7.82 NLOS
elseif nu_real_d2d == 3.19
    sigma_d2d = 7.82; % 4dB for LOS, 7.82 NLOS
end

%% USER real position
index_counter = 1;
for index_run = 1:runs
%% add the cycle
    for index_user = 1:number_of_evaluation_points
        r2 = Rd*sqrt(rand(1,1));
        theta2 = (pi/2)*rand(1,1);
        realX = (r2.*cos(theta2))'; % x coordinate of the known position
        realY = (r2.*sin(theta2))'; %y coordinate of the known position
        realX_array(1,index_user) = realX;
        realY_array(1,index_user) = realY;
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
                [estimanted_pos_cel(index_user,:)] = estimated_coordinate_RSS_cellular(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,Blockage_real_cel,Blockage_est_cel,nu_real_cel,nu_est_cel,sigma_cel);         
                error_cel(index_counter,:) = real(sqrt((estimanted_pos_cel(index_user,1) - realX).^2 + (estimanted_pos_cel(index_user,2) - realY).^2)); %% 

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
                
                GDoP_cellular = (inv(H'*H)); % equation 17
                GDoP_cellular_xy(index_counter,:) = [GDoP_cellular(1,1) GDoP_cellular(2,2)];
                GDoP_cellular_value(index_counter,:) = sqrt(GDoP_cellular(1,1) + GDoP_cellular(2,2));

%               % equation 30 geometric lower bound:
                A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
                     -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
                     -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6

                b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                    pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                    pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];

                nu1 = normrnd(0,sigma_cel);
                nu2 = normrnd(0,sigma_cel);
                nu3 = normrnd(0,sigma_cel);
                nu4 = normrnd(0,sigma_cel);
                epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
                           nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
                           nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];

                b_bar = b + epsilon; % equation 5
                eq_8 = [1+(sigma_cel*pi_i2^2)/(sigma_cel*pi_i1^2) 1 1
                        1 1+(sigma_cel*pi_i3^2)/(sigma_cel*pi_i1^2) 1
                        1 1 1+(sigma_cel*pi_i4^2)/(sigma_cel*pi_i1^2)];
                sum_epsilon = 4*sigma_cel*pi_i1^2*eq_8;
%               sum_epsilon = eq_8;
                p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
                   
                pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x1 = (p_roof(1) - X_a1/pi_i1);
                lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
                    
                pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x2 = (p_roof(1) - X_a2/pi_i2);
                lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
                    
                pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x3 = (p_roof(1) - X_a3/pi_i3);
                lambda_y3 = (p_roof(2) - Y_a3/pi_i3);
    
                pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x4 = (p_roof(1) - X_a4/pi_i4);
                lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
                    
                H_bar = [lambda_x1 lambda_y1; 
                        lambda_x2 lambda_y2;
                        lambda_x3 lambda_y3;
                        lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                    
                GLB_cel = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% %                     GLB = sigma_cel*inv(H_bar*H_bar'); %equation 29
                GLB_cel_xy(index_counter,:) = [GLB_cel(1,1) GLB_cel(2,2)];
                GLB_cel_value(index_counter,:) = sqrt(GLB_cel(1,1) + GLB_cel(2,2));                
                
            elseif mode == 2 %RSS devices
                % user distribution in a circle of radius Rd
                if mode_d2d == 1
                %% random in all the area

                    r2 = Rd*sqrt(rand(number_of_anchors,1));
                    theta2 = (pi/2)*rand(number_of_anchors,1);
                    X_B = (r2.*cos(theta2))';
                    Y_B = (r2.*sin(theta2))';    
    % %               plot(X_B,Y_B,'o');
                    
                    X = [X_B' Y_B']; 
                    X_a1 = X_B(1);
                    X_a2 = X_B(2);
                    X_a3 = X_B(3);
                    X_a4 = X_B(4);
                    Y_a1 = Y_B(1);
                    Y_a2 = Y_B(2);
                    Y_a3 = Y_B(3);
                    Y_a4 = Y_B(4);
% %                      

                    realXB_array(index_user,:) = X_B;
                    realYB_array(index_user,:) = Y_B;

                elseif mode_d2d==2

                    %% Random AROUND USER    
                    theta2 = rand(number_of_anchors,1)*2*pi;
                    X_B = realX + d0.*cos(theta2);
                    Y_B = realY + d0.*sin(theta2);

%                      X = [X_B' Y_B']; 
                    X_a1 = X_B(1);
                    X_a2 = X_B(2);
                    X_a3 = X_B(3);
                    X_a4 = X_B(4);
                    Y_a1 = Y_B(1);
                    Y_a2 = Y_B(2);
                    Y_a3 = Y_B(3);
                    Y_a4 = Y_B(4);
                    
                    realXB_array(index_user,:) = X_B;
                    realYB_array(index_user,:) = Y_B;
                    X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];

                elseif mode_d2d == 3
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
 
                elseif mode_d2d == 4
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

                [estimanted_pos_D2D(index_user,:)] = estimated_coordinate_RSS_D2D(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,Blockage_real_d2d,Blockage_est_d2d,nu_real_d2d,nu_est_d2d,sigma_d2d);
                error_D2D(index_counter,:) = real(sqrt((estimanted_pos_D2D(index_user,1) - realX).^2 + (estimanted_pos_D2D(index_user,2) - realY).^2)); %% 

                [estimanted_pos_D2D_HD(index_user,:)] = estimated_coordinate_RSS_D2D_HD(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,Blockage_real_d2d,Blockage_est_d2d,nu_real_d2d,nu_est_d2d,sigma_d2d);
                error_D2D_HD(index_counter,:) = real(sqrt((estimanted_pos_D2D_HD(index_user,1) - realX).^2 + (estimanted_pos_D2D_HD(index_user,2) - realY).^2)); %% 

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
                    
                GDoP_d2d = (inv(H'*H)); % equation 17
                GDoP_d2d_xy(index_counter,:) = [GDoP_d2d(1,1) GDoP_d2d(2,2)];
                GDoP_d2d_value(index_counter,:) = sqrt(GDoP_d2d(1,1) + GDoP_d2d(2,2));
% 
%                    % equation 30 geometric lower bound:
                A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
                     -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
                     -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6

                b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                     pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                     pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];

                nu1 = normrnd(0,sigma_d2d);
                nu2 = normrnd(0,sigma_d2d);
                nu3 = normrnd(0,sigma_d2d);
                nu4 = normrnd(0,sigma_d2d);
                epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
                           nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
                           nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];

                b_bar = b + epsilon; % equation 5
                eq_8 = [1+(sigma_d2d*pi_i2^2)/(sigma_d2d*pi_i1^2) 1 1
                        1 1+(sigma_d2d*pi_i3^2)/(sigma_d2d*pi_i1^2) 1
                        1 1 1+(sigma_d2d*pi_i4^2)/(sigma_d2d*pi_i1^2)];
                sum_epsilon = 4*sigma_d2d*pi_i1^2*eq_8;
%                    sum_epsilon = eq_8;
                p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
               
                pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x1 = (p_roof(1) - X_a1/pi_i1);
                lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
                    
                pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x2 = (p_roof(1) - X_a2/pi_i2);
                lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
                
                pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x3 = (p_roof(1) - X_a3/pi_i3);
                lambda_y3 = (p_roof(2) - Y_a3/pi_i3);

                pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x4 = (p_roof(1) - X_a4/pi_i4);
                lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
                
                H_bar = [lambda_x1 lambda_y1; 
                    lambda_x2 lambda_y2;
                    lambda_x3 lambda_y3;
                    lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                
                GLB_d2d = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% %                     GLB = sigma_cel*inv(H_bar*H_bar'); %equation 29
                GLB_d2d_xy(index_counter,:) = [GLB_d2d(1,1) GLB_d2d(2,2)];
                GLB_d2d_value(index_counter,:) = sqrt(GLB_d2d(1,1) + GLB_d2d(2,2));

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

                [estimanted_pos_RIS(index_user,:)] = estimated_coordinate_RSS_RIS(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN,Blockage_real_ris,Blockage_est_ris,nu_real_ris,nu_est_ris,sigma_ris);
                error_RIS(index_counter,:) = real(sqrt((estimanted_pos_RIS(index_user,1) - realX).^2 + (estimanted_pos_RIS(index_user,2) - realY).^2)); %% 
                
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
                    
                GDoP_RIS = (inv(H'*H)); % equation 17
                GDoP_RIS_xy(index_counter,:) = [GDoP_RIS(1,1) GDoP_RIS(2,2)];
                GDoP_RIS_value(index_counter,:) = sqrt(GDoP_RIS(1,1) + GDoP_RIS(2,2));

%                    % equation 30 geometric lower bound:
                A = [-2*(X_a2-X_a1) -2*(Y_a2-Y_a1);
                     -2*(X_a3-X_a1) -2*(Y_a3-Y_a1);
                     -2*(X_a4-X_a1) -2*(Y_a4-Y_a1)]; % equation 6

                b = [pi_i2^2 - X_a2^2 - Y_a2^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                     pi_i3^2 - X_a3^2 - Y_a3^2 - pi_i1^2 + X_a1^2 + Y_a1^2;
                     pi_i4^2 - X_a4^2 - Y_a4^2 - pi_i1^2 + X_a1^2 + Y_a1^2];

                nu1 = normrnd(0,sigma_ris);
                nu2 = normrnd(0,sigma_ris);
                nu3 = normrnd(0,sigma_ris);
                nu4 = normrnd(0,sigma_ris);
                epsilon = [nu2^2 + 2*pi_i2*nu2 - nu1^2 + 2*pi_i1*nu1;
                           nu3^2 + 2*pi_i3*nu3 - nu1^2 + 2*pi_i1*nu1;
                           nu4^2 + 2*pi_i4*nu4 - nu1^2 + 2*pi_i1*nu1];

                b_bar = b + epsilon; % equation 5
                eq_8 = [1+(sigma_ris*pi_i2^2)/(sigma_ris*pi_i1^2) 1 1
                        1 1+(sigma_ris*pi_i3^2)/(sigma_ris*pi_i1^2) 1
                        1 1 1+(sigma_ris*pi_i4^2)/(sigma_ris*pi_i1^2)];
                sum_epsilon = 4*sigma_ris*pi_i1^2*eq_8;
%                    sum_epsilon = eq_8;
                p_roof = inv(A'*inv(sum_epsilon)*A)*A'*inv(sum_epsilon)*b_bar; % equation 9
                   
                pi_i1 = real(sqrt((X_a1- p_roof(1))^2 + (Y_a1- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x1 = (p_roof(1) - X_a1/pi_i1);
                lambda_y1 = (p_roof(2) - Y_a1/pi_i1);
                    
                pi_i2 = real(sqrt((X_a2- p_roof(1))^2 + (Y_a2- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x2 = (p_roof(1) - X_a2/pi_i2);
                lambda_y2 = (p_roof(2) - Y_a2/pi_i2);
                
                pi_i3 = real(sqrt((X_a3- p_roof(1))^2 + (Y_a3- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x3 = (p_roof(1) - X_a3/pi_i3);
                lambda_y3 = (p_roof(2) - Y_a3/pi_i3);

                pi_i4 = real(sqrt((X_a4- p_roof(1))^2 + (Y_a4- p_roof(2))^2)); %real distance between target node and anchor
                lambda_x4 = (p_roof(1) - X_a4/pi_i4);
                lambda_y4 = (p_roof(2) - Y_a4/pi_i4);
                    
                H_bar = [lambda_x1 lambda_y1; 
                    lambda_x2 lambda_y2;
                    lambda_x3 lambda_y3;
                    lambda_x4 lambda_y4]; %Jacobian matrix of the measurements
                
                GLB_RIS = inv(H_bar'*inv(covariance_matrix)*H_bar);% equation 30
% %                     GLB = sigma_ris*inv(H_bar*H_bar'); %equation 29
                GLB_RIS_xy(index_counter,:) = [GLB_RIS(1,1) GLB_RIS(2,2)];
                GLB_RIS_value(index_counter,:) = sqrt(GLB_RIS(1,1) + GLB_RIS(2,2));                 
            
            end
        end
        
        %fuse estimates
        estimanted_pos_cel_D2D(index_user,:) = [ mean([estimanted_pos_cel(index_user,1) estimanted_pos_D2D(index_user,1)]), mean([estimanted_pos_cel(index_user,2) estimanted_pos_D2D(index_user,2)])];
        error_cel_D2D(index_counter,:) = real(sqrt((estimanted_pos_cel_D2D(index_user,1) - realX).^2 + (estimanted_pos_cel_D2D(index_user,2) - realY).^2)); %% 

        estimanted_pos_cel_D2D_HD(index_user,:) = [ mean([estimanted_pos_cel(index_user,1) estimanted_pos_D2D_HD(index_user,1)]), mean([estimanted_pos_cel(index_user,2) estimanted_pos_D2D_HD(index_user,2)])];
        error_cel_D2D_HD(index_counter,:) = real(sqrt((estimanted_pos_cel_D2D_HD(index_user,1) - realX).^2 + (estimanted_pos_cel_D2D_HD(index_user,2) - realY).^2)); %% 

        estimanted_pos_cel_RIS(index_user,:) = [ mean([estimanted_pos_cel(index_user,1) estimanted_pos_RIS(index_user,1)]), mean([estimanted_pos_cel(index_user,2) estimanted_pos_RIS(index_user,2)])];
        error_cel_RIS(index_counter,:) = real(sqrt((estimanted_pos_cel_RIS(index_user,1) - realX).^2 + (estimanted_pos_cel_RIS(index_user,2) - realY).^2)); %% 

        estimanted_pos_cel_D2D_RIS(index_user,:) = [ mean([estimanted_pos_cel(index_user,1) estimanted_pos_D2D(index_user,1) estimanted_pos_RIS(index_user,1)]), mean([estimanted_pos_cel(index_user,2) estimanted_pos_D2D(index_user,2) estimanted_pos_RIS(index_user,2)])];
        error_cel_D2D_RIS(index_counter,:) = real(sqrt((estimanted_pos_cel_D2D_RIS(index_user,1) - realX).^2 + (estimanted_pos_cel_D2D_RIS(index_user,2) - realY).^2)); %% 

        estimanted_pos_cel_D2D_RIS_HD(index_user,:) = [ mean([estimanted_pos_cel(index_user,1) estimanted_pos_D2D_HD(index_user,1) estimanted_pos_RIS(index_user,1)]), mean([estimanted_pos_cel(index_user,2) estimanted_pos_D2D_HD(index_user,2) estimanted_pos_RIS(index_user,2)])];
        error_cel_D2D_RIS_HD(index_counter,:) = real(sqrt((estimanted_pos_cel_D2D_RIS_HD(index_user,1) - realX).^2 + (estimanted_pos_cel_D2D_RIS_HD(index_user,2) - realY).^2)); %% 

        index_counter=index_counter+1;
    end

end      
%% CDF of INDIVIDUAL ERROR
%% stopped here
figure(2)
[h_cel,stats_cel] = cdfplot(error_cel);
hold on
[h_d2d,stats_d2d] = cdfplot(error_D2D);
hold on
[h_d2d_HD,stats_d2d_HD] = cdfplot(error_D2D_HD);
hold on
[h_RIS, stats_RIS] = cdfplot(error_RIS);
hold on
[h_cel_d2d,stats_cel_d2d] = cdfplot(error_cel_D2D);
hold on
[h_cel_d2d_HD,stats_cel_d2d_HD] = cdfplot(error_cel_D2D_HD);
hold on
[h_cel_RIS,stats_cel_RIS] = cdfplot(error_cel_RIS);
hold on
[h_cel_d2d_RIS,stats_cel_d2d_RIS] = cdfplot(error_cel_D2D_RIS);
hold on
[h_cel_d2d_RIS_HD,stats_cel_d2d_RIS_HD] = cdfplot(error_cel_D2D_RIS_HD);

srt_cel = {strcat('Cellular, \mu=',num2str(stats_cel.mean,'%.2f'),', \sigma=',num2str(stats_cel.std,'%.2f'))};
srt_d2d = {strcat('D2D FD, \mu=',num2str(stats_d2d.mean,'%.2f'), ', \sigma=',num2str(stats_d2d.std,'%.2f'))};
srt_d2d_HD = {strcat('D2D HD, \mu=',num2str(stats_d2d_HD.mean,'%.2f'), ', \sigma=',num2str(stats_d2d_HD.std,'%.2f'))};
srt_RIS = {strcat('RIS, \mu=',num2str(stats_RIS.mean,'%.2f'),', \sigma=',num2str(stats_RIS.std,'%.2f'))};
srt_cel_d2d = {strcat('Cellular+D2D FD, \mu=',num2str(stats_cel_d2d.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d.std,'%.2f'))};
srt_cel_d2d_HD = {strcat('Cellular+D2D HD, \mu=',num2str(stats_cel_d2d_HD.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_HD.std,'%.2f'))};
srt_cel_RIS = {strcat('Cellular+RIS, \mu=',num2str(stats_cel_RIS.mean,'%.2f'),', \sigma=',num2str(stats_cel_RIS.std,'%.2f'))};
srt_cel_d2d_RIS = {strcat('Cellular+D2D FD+RIS, \mu=',num2str(stats_cel_d2d_RIS.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_RIS.std,'%.2f'))};
srt_cel_d2d_RIS_HD = {strcat('Cellular+D2D HD+RIS, \mu=',num2str(stats_cel_d2d_RIS_HD.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_RIS_HD.std,'%.2f'))};
legend(srt_cel{1,1},srt_d2d{1,1},srt_d2d_HD{1,1},srt_RIS{1,1},srt_cel_d2d{1,1},srt_cel_d2d_HD{1,1},srt_cel_RIS{1,1},srt_cel_d2d_RIS{1,1},srt_cel_d2d_RIS_HD{1,1},'Location','best')
ylabel('CDF of individual error')
xlabel('Error, m')
set(h_cel,'LineWidth',1.5)
set(h_d2d,'LineWidth',1.5)
set(h_RIS,'LineWidth',1.5)
set(h_cel_d2d,'LineWidth',1.5)
set(h_cel_RIS,'LineWidth',1.5)
set(h_cel_d2d_RIS,'LineWidth',1.5)
grid on
hold off
xlim([0,600])

% 50th 75th 95th percentile
error_cel_prc = prctile(error_cel,[50 75 95]);
error_D2D_prc = prctile(error_D2D,[50 75 95]);
error_D2D_prc_HD = prctile(error_D2D_HD,[50 75 95]);
error_RIS_prc = prctile(error_RIS,[50 75 95]);
error_cel_D2D_prc = prctile(error_cel_D2D,[50 75 95]);
error_cel_D2D_prc_HD = prctile(error_cel_D2D_HD,[50 75 95]);
error_cel_RIS_prc = prctile(error_cel_RIS,[50 75 95]);
error_cel_D2D_RIS_prc = prctile(error_cel_D2D_RIS,[50 75 95]);
error_cel_D2D_RIS_prc_HD = prctile(error_cel_D2D_RIS_HD,[50 75 95]);


[h_CRLB_cellular,stats_CRLB_cellular] = cdfplot(CRLB_cellular_value);
hold on
[h_CRLB_d2d,stats_CRLB_d2d] = cdfplot(CRLB_d2d_value);
hold on
[h_CRLB_RIS,stats_CRLB_RIS] = cdfplot(CRLB_RIS_value);
hold on
[h_GDoP_cellular, stats_GDoP_cellular] = cdfplot(GDoP_cellular_value);
hold on
[h_GDoP_d2d,stats_GDoP_d2d] = cdfplot(GDoP_d2d_value);
hold on
[h_GDoP_RIS,stats_GDoP_RIS] = cdfplot(GDoP_RIS_value);
hold on
[h_GLB_cel,stats_GLB_cel] = cdfplot(GLB_cel_value);
hold on
[h_GLB_d2d,stats_GLB_d2d] = cdfplot(GLB_d2d_value);
hold on
[h_GLB_RIS,stats_GLB_RIS] = cdfplot(GLB_RIS_value);

% lower bounds
error_CRLB_cel_prc = prctile(CRLB_cellular_value ,[50 75 95]);
error_CRLB_d2d_prc = prctile(CRLB_d2d_value ,[50 75 95]);
error_CRLB_RIS_prc = prctile(CRLB_RIS_value ,[50 75 95]);
error_GDoP_cel_prc = prctile(GDoP_cellular_value ,[50 75 95]);
error_GDoP_d2d_prc = prctile(GDoP_d2d_value ,[50 75 95]);
error_GDoP_RIS_prc = prctile(GDoP_RIS_value ,[50 75 95]);
error_GLB_cel_prc = prctile(GLB_cel_value ,[50 75 95]);
error_GLB_d2d_prc = prctile(GLB_d2d_value ,[50 75 95]);
error_GLB_RIS_prc = prctile(GLB_RIS_value ,[50 75 95]);


%MSE RMSE
error_cel_MSE = mean(error_cel.^2);
error_D2D_MSE = mean(error_D2D.^2);
error_D2D_MSE_HD = mean(error_D2D_HD.^2);
error_RIS_MSE = mean(error_RIS.^2);
error_cel_D2D_MSE = mean(error_cel_D2D.^2);
error_cel_D2D_MSE_HD = mean(error_cel_D2D_HD.^2);
error_cel_RIS_MSE = mean(error_cel_RIS.^2);
error_cel_D2D_RIS_MSE = mean(error_cel_D2D_RIS.^2);
error_cel_D2D_RIS_MSE_HD = mean(error_cel_D2D_RIS_HD.^2);

error_CRLB_cel_MSE = mean(CRLB_cellular_value.^2);
error_CRLB_d2d_MSE = mean(CRLB_d2d_value.^2);
error_CRLB_RIS_MSE = mean(CRLB_RIS_value.^2);
error_GDoP_cel_MSE = mean(GDoP_cellular_value.^2);
error_GDoP_d2d_MSE = mean(GDoP_d2d_value.^2);
error_GDoP_RIS_MSE = mean(GDoP_RIS_value.^2);
error_GLB_cel_MSE = mean(GLB_cel_value.^2);
error_GLB_d2d_MSE = mean(GLB_d2d_value.^2);
error_GLB_RIS_MSE = mean(GLB_RIS_value.^2);

%RMSE
error_cel_RMSE = sqrt(error_cel_MSE);
error_D2D_RMSE = sqrt(error_D2D_MSE);
error_D2D_RMSE_HD = sqrt(error_D2D_MSE_HD);
error_RIS_RMSE = sqrt(error_RIS_MSE);
error_cel_D2D_RMSE = sqrt(error_cel_D2D_MSE);
error_cel_D2D_RMSE_HD = sqrt(error_cel_D2D_MSE_HD);
error_cel_RIS_RMSE = sqrt(error_cel_RIS_MSE);
error_cel_D2D_RIS_RMSE = sqrt(error_cel_D2D_RIS_MSE);
error_cel_D2D_RIS_RMSE_HD = sqrt(error_cel_D2D_RIS_MSE_HD);

error_CRLB_cel_RMSE = sqrt(error_CRLB_cel_MSE);
error_CRLB_d2d_RMSE = sqrt(error_CRLB_d2d_MSE);
error_CRLB_RIS_RMSE = sqrt(error_CRLB_RIS_MSE);
error_GDoP_cel_RMSE = sqrt(error_GDoP_cel_MSE);
error_GDoP_d2d_RMSE = sqrt(error_GDoP_d2d_MSE);
error_GDoP_RIS_RMSE = sqrt(error_GDoP_RIS_MSE);
error_GLB_cel_RMSE = sqrt(error_GLB_cel_MSE);
error_GLB_d2d_RMSE = sqrt(error_GLB_d2d_MSE);
error_GLB_RIS_RMSE = sqrt(error_GLB_RIS_MSE);

% save('bounds_d2d_random_100_by_100.mat');
% save('bounds_d2d_equal_placement.mat');
% save('bounds_d2d_20m.mat');
