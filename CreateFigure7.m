%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code creates Figure 7
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 2.0 (Last edited: 2023-03-16)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  load('RIS_all_LOS.mat'); 

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
 set(h_RIS1,'LineWidth',1.5, 'Color',[31/255 119/255 180/255])
 set(h_RIS2,'LineWidth',1.5, 'Color',[255/255 127/255 14/255 ])
 set(h_RIS3,'LineWidth',1.5, 'Color',[44/255 160/255 44/255 ])
 set(h_RIS4,'LineWidth',1.5, 'Color',[214/255 39/255 40/255])
 set(h_RIS5,'LineWidth',1.5, 'Color',[148/255 103/255 189/255])
 set(h_RIS6,'LineWidth',1.5,'Color',[0 204/255 204/255])
 set(h_RIS7,'LineWidth',1.5, 'Color',[96/255 96/255 96/255])
 set(h_RIS8,'LineWidth',1.5,'Color',[0/255 153/255 153/255])
 set(h_RIS9,'LineWidth',1.5, 'Color',[190/255 190/255 190/255])



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
 set(h_snr1,'LineWidth',1.5, 'Color',[31/255 119/255 180/255])
 set(h_snr2,'LineWidth',1.5, 'Color',[255/255 127/255 14/255 ])
 set(h_snr3,'LineWidth',1.5, 'Color',[44/255 160/255 44/255 ])
 set(h_snr4,'LineWidth',1.5, 'Color',[214/255 39/255 40/255])
 set(h_snr5,'LineWidth',1.5, 'Color',[148/255 103/255 189/255])
 set(h_snr6,'LineWidth',1.5,'Color',[0 204/255 204/255])
 set(h_snr7,'LineWidth',1.5, 'Color',[96/255 96/255 96/255])
 set(h_snr8,'LineWidth',1.5,'Color',[0/255 153/255 153/255])
 set(h_snr9,'LineWidth',1.5, 'Color',[190/255 190/255 190/255])
 grid on
 hold off

