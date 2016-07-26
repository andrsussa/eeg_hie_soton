clear all
close all

subjects = {'2.1.1', '2.4.2', '2.5.1', '2.5.2', '2.2.1', '2.2.2', '2.3.2'};

N = size(subjects,2);
m = 900;

% for i = 1:size(subjects,2) 
%     sub = cell2mat(subjects(i));
%     dotlessSub = ['sub' regexprep(sub, {'\.'},{''})];
%     EEG = pop_biosig(['H:\repos\eeg_hie_soton\RAWEEGDATA\' sub ' 1min']);
%     if ispc
% %     EEG = pop_biosig('H:\repos\eeg_hie_soton\2.1.1min');
%         EEG=pop_chanedit(EEG, 'lookup',['H:\\thesis\\eeglab13_5_4b\\plugins'...
%             '\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp']);
%     else
% %     EEG = pop_biosig('/home/andres/repos/eeg_hie_soton/2.1.1min');
%         EEG = pop_chanedit(EEG, 'lookup',['/home/andres/MATLAB/'...
%             'eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/'...
%             'standard-10-5-cap385.elp']);
%     end
%     EEG = pop_select(EEG, 'nochannel', {'FPZ' 'AUX1' 'AUX2' 'AUX3' 'AUX4',...
%         'AUX5' 'AUX6' 'AUX7' 'AUX8' 'PG1' 'PG2' 'A1' 'A2'});
%     eegplot(EEG.data, 'title', sub);
% end

epochs = zeros(1, size(subjects,2));
if exist('features.mat', 'file')
    load features.mat
else
    for i = 1:size(subjects,2)
        sub = cell2mat(subjects(i));
        dotlessSub = ['sub' regexprep(sub, {'\.'},{''})];
        [features.(dotlessSub).data, features.(dotlessSub).index, epochs(i)] =...
            eeg_hie_feaX(...
            ['H:\repos\eeg_hie_soton\RAWEEGDATA\' sub ' 1min']);
    end
    display('All subjects were successfully analysed.');
    name = 'features.mat';
    save(name,'features');
end

X = zeros(m,N);
for i = 1:N
    sub = cell2mat(subjects(i));
    dotlessSub = ['sub' regexprep(sub, {'\.'},{''})];
    X(:,i) = features.(dotlessSub).data;
end

N1 = 1:4;
N2 = 5:7;
X1 = zscore(X(:,N1), [], 1);
X2 = zscore(X(:,N2), [], 1);

FDR = ((mean(X1,2) - mean(X2,2)).^2) ./...
    ((var(X1,0,2).^2) + (var(X2,0,2).^2));

indexFeatures = features.sub211.index';
[~, inFDRALL] = sort(FDR,'descend');

topFEA = inFDRALL(1:7);
[hE5, pE] = ttest2(X1, X2, 'Vartype', 'equal', 'Dim', 2);
[hE1, ~] = ttest2(X1, X2, 'Vartype', 'equal', 'Dim', 2, 'Alpha', 0.01);
[hU5, pU] = ttest2(X1, X2, 'Vartype', 'unequal', 'Dim', 2);
[hU1, ~] = ttest2(X1, X2, 'Vartype', 'unequal', 'Dim', 2, 'Alpha', 0.01);
pWW = zeros(900,1);
for i = 1:900
    pWW(i) = ranksum(X1(i,:),X2(i,:));
end
hWW5 = pWW < 0.05;
hWW1 = pWW < 0.01;
sortedFeatures = table(indexFeatures(topFEA), FDR(topFEA), pE(topFEA),...
    hE5(topFEA), hE1(topFEA), pU(topFEA), hU5(topFEA), hU1(topFEA),...
    pWW(topFEA), hWW5(topFEA), hWW1(topFEA), 'VariableNames',...
    {'Index' 'FDR' 'pvalueEV' 'hypoEV5' 'hypoEV1' 'pvalueUV' 'hypoUV5'...
    'hypoUV1' 'pvalueWW' 'hypoWW5' 'hypoWW1'});

X1t = X1(topFEA,:);
X2t = X2(topFEA,:);
Xt = [X1t';X2t'];
Y=[1*ones(1,4) 2*ones(1,3)];

figure(1);
for i = 1:size(Xt,2)
    Xi = Xt(:,i);
    MdlLinear = fitcdiscr(Xi,Y');

    L = MdlLinear.Coeffs(1,2).Linear;
    t1=L'*Xi(Y==1,:)';
    t2=L'*Xi(Y==2,:)';

    figure(1), subplot(2,4,i), h1 = histfit(t1);
    figure(1), hold on
    figure(1), subplot(2,4,i), h2 = histfit(t2);
    h1(1).FaceColor = [1 .8 .8];
    h2(1).FaceColor = [.8 .8 1];
    h2(2).Color = [0 0 1];
end

corrplot(Xt);

% Xt(:,2) = [];
% Xt(:,6-1) = [];
Xt(:,6:7) = [];

MdlLinear = fitcdiscr(Xt,Y');

L = MdlLinear.Coeffs(1,2).Linear;
t1=L'*Xt(Y==1,:)';
t2=L'*Xt(Y==2,:)';

figure(3), h1 = histfit(t1);
figure(3), hold on
figure(3), h2 = histfit(t2);
h1(1).FaceColor = [1 .8 .8];
h2(1).FaceColor = [.8 .8 1];
h2(2).Color = [0 0 1];