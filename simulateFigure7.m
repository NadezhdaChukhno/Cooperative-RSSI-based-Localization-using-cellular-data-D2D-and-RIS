%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig.7 CDF of individual error, SNR for RIS with% different %
% number of antenna elements, LOS, PLE=2.1, SF=4.                                %                              %
% Article: [D2D-aided versus RIS-aided Cooperative Positioning: Theoretical Model% 
% for RSSI-based Ranging and Performance Comparison]                             % 
% Download article: [link]                                                       %
% This is version 1.0 (Last edited: 2022-07-14)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all 
%% RSS  location estimation
%% parameters
%% BS/RIS/USER/D2D coordinates
% 4 coordinates of known anchors
number_of_anchors=4;
Rd=100; % 100 x 100 square
mode_array=[3];
% mode=3;
%1 - cellular
%2 - D2D
%3 - RIS
NN=1024; % Number of reflective elements
mode_RIS_array=[1,2,3,4,5,6,7,8,9];

%% USER real position
number_of_evaluation_points=30;
runs=1000;
index_counter=1;
for index_run=1:runs
%% add the cycle
    for index_user=1:number_of_evaluation_points
        r2 = Rd*sqrt(rand(1,1));
        theta2 = (pi/2)*rand(1,1);
        realX= (r2.*cos(theta2))'; % x coordinate of the known position
        realY= (r2.*sin(theta2))'; %y coordinate of the known position
        realX_array(1,index_user)=realX;
        realY_array(1,index_user)=realY;
        %% 
        for index_mode = 1:length(mode_array)
            mode= mode_array(index_mode);
            if mode==3 %% RISs
                % RISs are located clode to the BSs
                X_a1=10;
                Y_a1=10;
                X_a2=10;
                Y_a2=90;
                X_a3=90;
                Y_a3=10;
                X_a4=90;
                Y_a4=90;
                X = [X_a1, Y_a1; X_a2, Y_a2; X_a3,Y_a3; X_a4,Y_a4];
            %     plot (X_a1,Y_a1,'o',X_a2,Y_a2,'*',X_a3,Y_a3,'s',X_a4,Y_a4,'^')

                %BS positions
                X_a1_BS=0;
                Y_a1_BS=0;
                X_a2_BS=0;
                Y_a2_BS=100;
                X_a3_BS=100;
                Y_a3_BS=0;
                X_a4_BS=100;
                Y_a4_BS=100;

                    mu=0;
                    sigma=4; % 4dB for LOS, 7,82 NLOS
                    w_dB_real = normrnd(mu,sigma); % noise for NLOS
                    w_dB_est = normrnd(mu,sigma);
                for index_RIS_mode = 1:length(mode_RIS_array)
                    mode_RIS= mode_RIS_array(index_RIS_mode);

                    if mode_RIS==1
                        NN=8;
                        [estimanted_pos_RIS1(index_user,:),snr_1(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS1(index_counter,:)=real(sqrt( (estimanted_pos_RIS1(index_user,1)- realX).^2 + ( estimanted_pos_RIS1(index_user,2)- realY).^2 )); %% 
                        snr1(index_counter,:)=snr_1(index_user,:);
                    elseif  mode_RIS==2
                        NN=16;
                        [estimanted_pos_RIS2(index_user,:),snr_2(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS2(index_counter,:)=real(sqrt( (estimanted_pos_RIS2(index_user,1)- realX).^2 + ( estimanted_pos_RIS2(index_user,2)- realY).^2 )); %% 
                        snr2(index_counter,:)=snr_2(index_user,:);
                    elseif  mode_RIS==3
                        NN=32;
                        [estimanted_pos_RIS3(index_user,:),snr_3(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS3(index_counter,:)=real(sqrt( (estimanted_pos_RIS3(index_user,1)- realX).^2 + ( estimanted_pos_RIS3(index_user,2)- realY).^2 )); %% 
                        snr3(index_counter,:)=snr_3(index_user,:);
                    elseif  mode_RIS==4
                        NN=64;
                        [estimanted_pos_RIS4(index_user,:),snr_4(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS4(index_counter,:)=real(sqrt( (estimanted_pos_RIS4(index_user,1)- realX).^2 + ( estimanted_pos_RIS4(index_user,2)- realY).^2 )); %% 
                        snr4(index_counter,:)=snr_4(index_user,:); 
                    elseif  mode_RIS==5
                        NN=128;
                        [estimanted_pos_RIS5(index_user,:),snr_5(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS5(index_counter,:)=real(sqrt( (estimanted_pos_RIS5(index_user,1)- realX).^2 + ( estimanted_pos_RIS5(index_user,2)- realY).^2 )); %% 
                        snr5(index_counter,:)=snr_5(index_user,:); 
                       
                    elseif  mode_RIS==6
                         NN=256;
                        [estimanted_pos_RIS6(index_user,:),snr_6(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS6(index_counter,:)=real(sqrt( (estimanted_pos_RIS6(index_user,1)- realX).^2 + ( estimanted_pos_RIS6(index_user,2)- realY).^2 )); %% 
                        snr6(index_counter,:)=snr_6(index_user,:);
                    elseif  mode_RIS==7
                        NN=512;
                        [estimanted_pos_RIS7(index_user,:),snr_7(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS7(index_counter,:)=real(sqrt( (estimanted_pos_RIS7(index_user,1)- realX).^2 + ( estimanted_pos_RIS7(index_user,2)- realY).^2 )); %% 
                        snr7(index_counter,:)=snr_7(index_user,:);
                    elseif  mode_RIS==8
                        NN=1024;
                        [estimanted_pos_RIS8(index_user,:),snr_8(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS8(index_counter,:)=real(sqrt( (estimanted_pos_RIS8(index_user,1)- realX).^2 + ( estimanted_pos_RIS8(index_user,2)- realY).^2 )); %% 
                        snr8(index_counter,:)=snr_8(index_user,:);
                    elseif  mode_RIS==9
                        NN=2048;
                        [estimanted_pos_RIS9(index_user,:),snr_9(index_user,:)]=estimated_coordinate_RSS_RIS_new(realX,realY,X_a1,Y_a1,X_a2,Y_a2,X_a3,Y_a3,X_a4,Y_a4,X,X_a1_BS,X_a2_BS,X_a3_BS,X_a4_BS,Y_a1_BS,Y_a2_BS,Y_a3_BS,Y_a4_BS,NN);
                        error_RIS9(index_counter,:)=real(sqrt( (estimanted_pos_RIS9(index_user,1)- realX).^2 + ( estimanted_pos_RIS9(index_user,2)- realY).^2 )); %% 
                        snr9(index_counter,:)=snr_9(index_user,:);    
                    end
                end
                
             end
        end
         index_counter=index_counter+1;
    end
 
end      

 figure
  [h_RIS1,stats_RIS1] = cdfplot(error_RIS1);
  hold on
  [h_RIS2,stats_RIS2] = cdfplot(error_RIS2);
  hold on
  [h_RIS3,stats_RIS3] = cdfplot(error_RIS3);
  hold on
  [h_RIS4,stats_RIS4] = cdfplot(error_RIS4);
  hold on
    [h_RIS5,stats_RIS5] = cdfplot(error_RIS5);
  hold on
    [h_RIS6,stats_RIS6] = cdfplot(error_RIS6);
  hold on
    [h_RIS7,stats_RIS7] = cdfplot(error_RIS7);
  hold on
    [h_RIS8,stats_RIS8] = cdfplot(error_RIS8);
  hold on
    [h_RIS9,stats_RIS9] = cdfplot(error_RIS9);
  hold on
 srt_d2d1={strcat('RIS 8, \mu=',num2str(stats_RIS1.mean,'%.2f'), ', \sigma=',num2str(stats_RIS1.std,'%.2f'))};
 srt_d2d2={strcat('RIS 16, \mu=',num2str(stats_RIS2.mean,'%.2f'), ', \sigma=',num2str(stats_RIS2.std,'%.2f'))};
 srt_d2d3={strcat('RIS 32, \mu=',num2str(stats_RIS3.mean,'%.2f'), ', \sigma=',num2str(stats_RIS3.std,'%.2f'))};
 srt_d2d4={strcat('RIS 64, \mu=',num2str(stats_RIS4.mean,'%.2f'), ', \sigma=',num2str(stats_RIS4.std,'%.2f'))};
srt_d2d5={strcat('RIS 128, \mu=',num2str(stats_RIS5.mean,'%.2f'), ', \sigma=',num2str(stats_RIS5.std,'%.2f'))};
 srt_d2d6={strcat('RIS 256, \mu=',num2str(stats_RIS6.mean,'%.2f'), ', \sigma=',num2str(stats_RIS6.std,'%.2f'))};
 srt_d2d7={strcat('RIS 512, \mu=',num2str(stats_RIS7.mean,'%.2f'), ', \sigma=',num2str(stats_RIS7.std,'%.2f'))};
 srt_d2d8={strcat('RIS 1024, \mu=',num2str(stats_RIS8.mean,'%.2f'), ', \sigma=',num2str(stats_RIS8.std,'%.2f'))};
 srt_d2d9={strcat('RIS 2048, \mu=',num2str(stats_RIS9.mean,'%.2f'), ', \sigma=',num2str(stats_RIS9.std,'%.2f'))};
 legend(srt_d2d1{1,1},srt_d2d2{1,1},srt_d2d3{1,1},srt_d2d4{1,1},srt_d2d5{1,1},srt_d2d6{1,1},srt_d2d7{1,1},srt_d2d8{1,1},srt_d2d9{1,1},'Location','best')
 ylabel('CDF of individual error')
 xlabel('Error, m')
 set(h_RIS1,'LineWidth',1.5)
 set(h_RIS2,'LineWidth',1.5)
 set(h_RIS3,'LineWidth',1.5)
 set(h_RIS4,'LineWidth',1.5)
 set(h_RIS5,'LineWidth',1.5)
 set(h_RIS6,'LineWidth',1.5)
 set(h_RIS7,'LineWidth',1.5)
 set(h_RIS8,'LineWidth',1.5)
 set(h_RIS9,'LineWidth',1.5)
 grid on
 xlim([0,600])
 hold off


 figure
  [h_snr1,stats_snr1] = cdfplot(snr1);
  hold on
  [h_snr2,stats_snr2] = cdfplot(snr2);
  hold on
  [h_snr3,stats_snr3] = cdfplot(snr3);
  hold on
  [h_snr4,stats_snr4] = cdfplot(snr4);
  hold on
   [h_snr5,stats_snr5] = cdfplot(snr5);
  hold on
  [h_snr6,stats_snr6] = cdfplot(snr6);
  hold on
  [h_snr7,stats_snr7] = cdfplot(snr7);
  hold on
  [h_snr8,stats_snr8] = cdfplot(snr8);
  hold on
  [h_snr9,stats_snr9] = cdfplot(snr9);
  hold on 
  srt_d2d1={strcat('RIS 8, \mu=',num2str(stats_snr1.mean,'%.2f'))};
 srt_d2d2={strcat('RIS 16, \mu=',num2str(stats_snr2.mean,'%.2f'))};
 srt_d2d3={strcat('RIS 32, \mu=',num2str(stats_snr3.mean,'%.2f'))};
 srt_d2d4={strcat('RIS 64, \mu=',num2str(stats_snr4.mean,'%.2f'))};
srt_d2d5={strcat('RIS 128, \mu=',num2str(stats_snr5.mean,'%.2f'))};
  srt_d2d6={strcat('RIS 256, \mu=',num2str(stats_snr6.mean,'%.2f'))};
 srt_d2d7={strcat('RIS 512, \mu=',num2str(stats_snr7.mean,'%.2f'))};
 srt_d2d8={strcat('RIS 1024, \mu=',num2str(stats_snr8.mean,'%.2f'))};
 srt_d2d9={strcat('RIS 2048, \mu=',num2str(stats_snr9.mean,'%.2f'))};
 legend(srt_d2d1{1,1},srt_d2d2{1,1},srt_d2d3{1,1},srt_d2d4{1,1},srt_d2d5{1,1},srt_d2d6{1,1},srt_d2d7{1,1},srt_d2d8{1,1},srt_d2d9{1,1},'Location','best')
 ylabel('CDF of SNR')
 xlabel('SNR, dBm')
 set(h_RIS1,'LineWidth',1.5)
 set(h_RIS2,'LineWidth',1.5)
 set(h_RIS3,'LineWidth',1.5)
 set(h_RIS4,'LineWidth',1.5)
 set(h_RIS5,'LineWidth',1.5)
 set(h_RIS6,'LineWidth',1.5)
 set(h_RIS7,'LineWidth',1.5)
 set(h_RIS8,'LineWidth',1.5)
 set(h_RIS9,'LineWidth',1.5)
 grid on
 hold off

%      save('EXP3_RIS_all_LOS2.mat'); 