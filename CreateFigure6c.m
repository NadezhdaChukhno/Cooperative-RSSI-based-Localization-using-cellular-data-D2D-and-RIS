%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code creates Figure 6c and Table Vc
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 2.0 (Last edited: 2023-03-16)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('figure6c.mat')
figure(1)
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

 srt_cel={strcat('Cellular, \mu=',num2str(stats_cel.mean,'%.2f'),', \sigma=',num2str(stats_cel.std,'%.2f'))};
 srt_d2d={strcat('D2D FD, \mu=',num2str(stats_d2d.mean,'%.2f'), ', \sigma=',num2str(stats_d2d.std,'%.2f'))};
 srt_d2d_HD={strcat('D2D HD, \mu=',num2str(stats_d2d_HD.mean,'%.2f'), ', \sigma=',num2str(stats_d2d_HD.std,'%.2f'))};
 srt_RIS={strcat('RIS, \mu=',num2str(stats_RIS.mean,'%.2f'),', \sigma=',num2str(stats_RIS.std,'%.2f'))};
 srt_cel_d2d={strcat('Cellular+D2D FD, \mu=',num2str(stats_cel_d2d.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d.std,'%.2f'))};
 srt_cel_d2d_HD={strcat('Cellular+D2D HD, \mu=',num2str(stats_cel_d2d_HD.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_HD.std,'%.2f'))};
 srt_cel_RIS={strcat('Cellular+RIS, \mu=',num2str(stats_cel_RIS.mean,'%.2f'),', \sigma=',num2str(stats_cel_RIS.std,'%.2f'))};
 srt_cel_d2d_RIS={strcat('Cellular+D2D FD+RIS, \mu=',num2str(stats_cel_d2d_RIS.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_RIS.std,'%.2f'))};
 srt_cel_d2d_RIS_HD={strcat('Cellular+D2D HD+RIS, \mu=',num2str(stats_cel_d2d_RIS_HD.mean,'%.2f'),', \sigma=',num2str(stats_cel_d2d_RIS_HD.std,'%.2f'))};
 legend(srt_cel{1,1},srt_d2d{1,1},srt_d2d_HD{1,1},srt_RIS{1,1},srt_cel_d2d{1,1},srt_cel_d2d_HD{1,1},srt_cel_RIS{1,1},srt_cel_d2d_RIS{1,1},srt_cel_d2d_RIS_HD{1,1},'Location','best')
 ylabel('CDF of individual error')
 xlabel('Error, m')
 set(h_cel,'LineWidth',1.5, 'Color',[31/255 119/255 180/255])
 set(h_d2d,'LineWidth',2.9, 'Color',[255/255 127/255 14/255])
 set(h_d2d_HD,'LineWidth',1.5, 'Color',[44/255 160/255 44/255])
 set(h_RIS,'LineWidth',1.5, 'Color',[214/255 39/255 40/255])
 set(h_cel_d2d,'LineWidth',2.9, 'Color',[148/255 103/255 189/255])
 set(h_cel_d2d_HD,'LineWidth',1.5,'Color',[0 204/255 204/255])
 set(h_cel_RIS,'LineWidth',1.5, 'Color',[96/255 96/255 96/255] )
 set(h_cel_d2d_RIS,'LineWidth',2.9,'Color',[0/255 153/255 153/255])
 set(h_cel_d2d_RIS_HD,'LineWidth',1.5, 'Color',[190/255 190/255 190/255])

 grid on
 hold off
 xlim([0,300])

 
% 50th 75th 95th percentile
error_cel_prc=prctile(error_cel,[50 75 95]);
error_D2D_prc=prctile(error_D2D,[50 75 95]);
error_D2D_prc_HD=prctile(error_D2D_HD,[50 75 95]);
error_RIS_prc=prctile(error_RIS,[50 75 95]);
error_cel_D2D_prc=prctile(error_cel_D2D,[50 75 95]);
error_cel_D2D_prc_HD=prctile(error_cel_D2D_HD,[50 75 95]);
error_cel_RIS_prc=prctile(error_cel_RIS,[50 75 95]);
error_cel_D2D_RIS_prc=prctile(error_cel_D2D_RIS,[50 75 95]);
error_cel_D2D_RIS_prc_HD=prctile(error_cel_D2D_RIS_HD,[50 75 95]);
%MSE RMSE

error_cel_MSE=mean(error_cel.^2);
error_D2D_MSE=mean(error_D2D.^2);
error_D2D_MSE_HD=mean(error_D2D_HD.^2);
error_RIS_MSE=mean(error_RIS.^2);
error_cel_D2D_MSE=mean(error_cel_D2D.^2);
error_cel_D2D_MSE_HD=mean(error_cel_D2D_HD.^2);
error_cel_RIS_MSE=mean(error_cel_RIS.^2);
error_cel_D2D_RIS_MSE=mean(error_cel_D2D_RIS.^2);
error_cel_D2D_RIS_MSE_HD=mean(error_cel_D2D_RIS_HD.^2);

error_cel_RMSE=sqrt(error_cel_MSE);
error_D2D_RMSE=sqrt(error_D2D_MSE);
error_D2D_RMSE_HD=sqrt(error_D2D_MSE_HD);
error_RIS_RMSE=sqrt(error_RIS_MSE);
error_cel_D2D_RMSE=sqrt(error_cel_D2D_MSE);
error_cel_D2D_RMSE_HD=sqrt(error_cel_D2D_MSE_HD);
error_cel_RIS_RMSE=sqrt(error_cel_RIS_MSE);
error_cel_D2D_RIS_RMSE=sqrt(error_cel_D2D_RIS_MSE);
error_cel_D2D_RIS_RMSE_HD=sqrt(error_cel_D2D_RIS_MSE_HD);