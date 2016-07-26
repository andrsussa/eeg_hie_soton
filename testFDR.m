clear all
close all

randn('seed',0) %For reproducibility of the results 
S=[.3 1.5; 1.5 9]; % Covariance matrix
[lx,ly]=size(S); 
mv=[-1 0; 1 0]'; % Mean for gaussian distribution
N=200; % Length of class
X2=[mvnrnd(mv(:,1),S,N); mvnrnd(mv(:,2),S,N)]'; 
y2=[ones(1,N), 2*ones(1,N)];

mv_est(:,1) = mean(X2(:,y2==1),2); 
mv_est(:,2) = mean(X2(:,y2==2),2);

[Sw,Sb,Sm]=scatter_mat(X2,y2);

w = Sw\(mv_est(:,1)-mv_est(:,2));

figure(1), plot(X2(1,y2==1),X2(2,y2==1),'r.',...
    X2(1,y2==2),X2(2,y2==2),'bo')
figure(1), axis equal 
%Computation of the projections 
t1=w'*X2(:,y2==1);
t2=w'*X2(:,y2==2);

figure(2), h1 = histfit(t1);
figure(2), hold on
figure(2), h2 = histfit(t2);
h1(1).FaceColor = [1 .8 .8];
h2(1).FaceColor = [.8 .8 1];
h2(2).Color = [0 0 1];

X_proj1=[t1;t1].*((w/(w'*w))*ones(1,length(t1))); 
X_proj2=[t2;t2].*((w/(w'*w))*ones(1,length(t2))); 

%Plot of the projections 
figure(1), hold on 
figure(1), plot(X_proj1(1,:),X_proj1(2,:),'y.',... 
    X_proj2(1,:),X_proj2(2,:),'co')