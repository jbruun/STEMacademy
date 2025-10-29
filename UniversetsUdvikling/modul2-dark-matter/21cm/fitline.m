clear all
clc
close all

filename = '21cm.data';
A = importdata(filename);
delimiterIn = ' ';
headerlinesIn = 3;

vel = A.data(:,1);
flux = A.data(:,2)+2.3;

figure
plot(vel,flux);

%Try to fit with the asymmetric Busy function from 
%https://arxiv.org/pdf/1311.5308.pdf
ft = fittype('a/4 * ( erf(b*(d+(x-e))) +1) * ( erf(c*(d-(x-e))) +1) * (f*(x-e)^2+1)', 'Independent', 'x');
startguess = [15.,0.001,0.01, 170.,2300.,0.015];
f = fit(vel,flux,ft, 'StartPoint', startguess);
cc = cfit(f);
conf = confint(f);
%Velocity
'Fitted central velocity from Simplified Busy function (km/s):';,cc.e
'95% Confidence interval from SBF (AA):',conf(:,1)
'95% Confidence interval from SBF (AA):',conf(:,2)
'95% Confidence interval from SBF (AA):',conf(:,3)
'95% Confidence interval from SBF (AA):',conf(:,4)
'95% Confidence interval from SBF (AA):',conf(:,5)
'95% Confidence interval from SBF (AA):',conf(:,6)

plot(f,vel,flux);
hold on
