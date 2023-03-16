%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig.1 LoS probability as a function of 2D distance between % 
% the BS and the MT according to 3GPP UMi Street Canyon model.                   %
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                  % 
% Download article: [link]                                                       %
% This is version 2.0 (Last edited: 2023-03-16)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
figure
d_array=0:5:600; % min distance:step:max distance
for i=1:length(d_array)
    d = d_array(i);
    if d<=18
        P_LOS_UMI(i)=1;
    else
        P_LOS_UMI(i)=18/d+exp(-d/36)*(1-18/d);
    end
end

plot(d_array,P_LOS_UMI*100,'LineWidth', 1.5, 'Color', [44/255 160/255 44/255])
ylabel('LoS probability, p_L, %')
xlabel('Distance between BS and MT, d_{2D}, m')
grid on