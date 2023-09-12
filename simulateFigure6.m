%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is code that generates plot with pdf of noise (Fig. 6).                                                        %
% Article: [Are D2D and RIS in the Same League? Cooperative RSSI-based 
% Localization Model and Performance Comparison]                                 % 
% Download article: [link]                                                       %
% This is version 3.0 (Last edited: 2023-09-06)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
figure (1)
mu = 0;
sigma = 4;
x = [-15:0.1:15];
y = normpdf(x,mu,sigma);
hFig = figure(1);
xwidth = 560;
ywidth = 260;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 xwidth ywidth])
plot(x,y,'LineWidth', 1.5, 'Color', [44/255 160/255 44/255])
xlabel('Possible noise values, dB')
ylabel('Probability density function')
str = {'\mu = 0'};
str2 = {'\sigma_{SF} = 4'};
text(10.7,0.095,str)
text(10.7,0.085,str2)
grid on