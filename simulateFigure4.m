%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code generates Fig.4 Visualisation of median and mean on skewed and       %
% normal distribution of data.                                                   %
% Article: [name]                                                                %
% Download article: [link]                                                       %
% This is version 1.0 (Last edited: 2022-07-13)                                  %
% Author: N. Chukhno                                                             %
% University Mediterranea of Reggio Calabria, Italy and CNIT, Italy.             %
% Universitat Jaume I, Spain                                                     %
% Email: nadezda.chukhno@unirc.it                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
%% skewed
figure
subplot(2,1,1)
rng default;  % For reproducibility
b = betarnd(2,10,100,1);
%Construct a histogram using 10 bins with a beta distribution fit.
histfit(b,10,'beta')
x1=xline(mean(b));
x1.LineWidth=2;
x1.Color='g';
x1.Label='mean';
x1.LabelOrientation='horizontal';
x1.LabelVerticalAlignment='middle';
x2=xline(median(b));
x2.LineWidth=2;
x2.Color='b';
x2.Label='median';
x2.LabelOrientation='horizontal';
x2.LabelHorizontalAlignment='left';

%% normal
subplot(2,1,2)
pd = makedist('Lognormal','mu',5,'sigma',2);
rng('default');  % For reproducibility
x = random(pd,100,1);
logx = log(x);
histfit(logx)
x1=xline(mean(logx));
x1.LineWidth=2;
x1.Color='g';
x1.Label='mean';
x1.LabelOrientation='horizontal';
x1.LabelVerticalAlignment='middle';
x2=xline(median(logx));
x2.LineWidth=2;
x2.Color='b';
x2.Label='median';
x2.LabelOrientation='horizontal';
x2.LabelHorizontalAlignment='left';

