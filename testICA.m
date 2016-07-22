clear all
close all

t = 0:1/100:100;

x1 = 2*sin(0.2*pi*t);
x2 = 2*sawtooth(0.1*pi*t,1);

x3 = [x1; x2];

plot(t, x1, t, x2);

[mixSingal,w,wi,rc] = instamix( [x1;x2], .01);

figure, plot(t, mixSingal);

[icasig, A, W] = fastica(mixSingal);

figure, plot(t, icasig);