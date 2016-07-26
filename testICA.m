clear all
close all

t = 0:1/100:100;

x1 = 2*sin(0.2*pi*t);
x2 = 2*sawtooth(0.1*pi*t,1);

fi = 1;
figure(fi), plot(t, x1, t, x2);

x3 = x1 - 2*x2;
x4 = 1.73*x1 +3.41*x2;

fi = fi + 1;
figure(fi), subplot(2,1,1), plot(t, x3);
figure(fi), subplot(2,1,2), plot(t, x4);

mixSignal = [x3; x4];
% [mixSignal,~,~,~] = instamix( [x1;x2], .01);

fi = fi + 1;
figure(fi), plot(t, mixSignal);

[weights,sphere] = runica(mixSignal);

runicaSig = weights \ mixSignal;

fi = fi + 1;
figure(fi), plot(t, runicaSig);

[icasig, A, W] = fastica(mixSignal);

fi = fi + 1;
figure(fi), plot(t, icasig);