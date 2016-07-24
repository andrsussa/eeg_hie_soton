function [ features, featuresIndex ] = eeg_hie_feaX( data_path )
%EEG_HIE_FEAX Extract features from EEG data.

EEG = pop_biosig(data_path);
if ispc
%     EEG = pop_biosig('H:\repos\eeg_hie_soton\2.1.1min');
    EEG=pop_chanedit(EEG, 'lookup',['H:\\thesis\\eeglab13_5_4b\\plugins'...
        '\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp']);
else
%     EEG = pop_biosig('/home/andres/repos/eeg_hie_soton/2.1.1min');
    EEG = pop_chanedit(EEG, 'lookup',['/home/andres/MATLAB/'...
        'eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/'...
        'standard-10-5-cap385.elp']);
end

%% Filtering
EEG = pop_select(EEG, 'nochannel', {'FPZ' 'AUX1' 'AUX2' 'AUX3' 'AUX4',...
    'AUX5' 'AUX6' 'AUX7' 'AUX8' 'PG1' 'PG2' 'A1' 'A2'});
notchlo = 49;
notchhi = 51;
bandlo = 0.5;
bandhi = 45;
epochlen = 5;
EEG = pop_eegfiltnew(EEG, notchlo, notchhi, [], 1, [], 0);
EEG = pop_eegfiltnew(EEG, bandlo, bandhi, [], 0, [], 0);

%% Epoching I
EEG.data(1,1:EEG.srate*epochlen:EEG.pnts) = 1;  % Set events
lim = [-1 2];
EEG = pop_chanevent(EEG, 1, 'edge', 'leading', 'edgelen', 0,...
    'duration', 'on');
EEG = pop_epoch(EEG, { }, lim, 'newname', '1MIN data epochs',...
    'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-1000 0]);   % Check what Baseline removing is 
                                    % for, and proper Value!!!!

%% Thresholding
chanNum = EEG.nbchan;
pnts = EEG.pnts;
thuV = 100;
EEG = pop_eegthresh(EEG, 1, 1:chanNum, -1*thuV, thuV, -1, 1.998, 0, 1);

%% ICA
EEGdata = reshape(EEG.data(1:chanNum,:,:), chanNum, EEG.pnts*EEG.trials);
EEGdata = double(EEGdata);
EEGica = fastica(EEGdata, 'stabilization', 'on');

%% Epoching II
srate = EEG.srate;
valuelim = [-Inf Inf];
tmpevent = EEG.event;
alllatencies = [ tmpevent(:).latency ];

EEGdata = epoch(EEGica, alllatencies, [lim(1) lim(2)]*srate,...
    'valuelim', valuelim, 'allevents', alllatencies);


%% Preparing structures
measuresCell = {'COH', 'iCOH', 'PLV', 'PLI', 'RHO'};
measuresDataCell = {[]; []; []; []; []};
measures = cell2struct(measuresDataCell, measuresCell, 1);
measuresNum = size(measuresCell,2);

bcMeasuresCell = {'Data', 'Modul', 'Trans', 'CharPath', 'Effi',...
    'NetRad', 'NetDia'};
bcMeasuresDataCell = {[]; []; []; []; []; []; []};
bcMeasures = cell2struct(bcMeasuresDataCell, bcMeasuresCell, 1);
% stabcMeasures = cell2struct(bcMeasuresDataCell(2:end), bcMeasuresCell(2:end), 1);
bcMeasuresNum = size(bcMeasuresCell,2) - 1;

staMeasuresCell = {'Mean', 'Median', 'StD', 'IQR',...
    'Skew', 'Kurt'};
% staMeasuresDataCell = {[]; []; []; []; []; []};
% staMeasures = cell2struct(staMeasuresDataCell, staMeasuresCell, 1);
staMeasuresNum = size(staMeasuresCell,2);

fBandsCell = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
bandsNum = size(fBandsCell,2);
for i = 1:bandsNum
    fBands.(cell2mat(fBandsCell(i))) = measures;
%     fFeatures.(cell2mat(fBandsCell(i))) = measures;
    for j = 1:measuresNum
        fBands.(cell2mat(fBandsCell(i))).(cell2mat(measuresCell(j))) = ...
            bcMeasures;
%         fFeatures.(cell2mat(fBandsCell(i))).(cell2mat(measuresCell(j)))...
%             = stabcMeasures;
%         for k = 2:bcMeasuresNum+1
%             fFeatures.(cell2mat(fBandsCell(i)))...
%                 .(cell2mat(measuresCell(j)))...
%                 .(cell2mat(bcMeasuresCell(k))) = staMeasures;
%         end
    end
        
end

clear measures measuresDataCell bcMeasuresDataCell bcMeasures

%% Window Properties
% winLength = EEG.pnts*1000/EEG.srate; % Lenght in ms
winLength = 1000;
baseline = 0;
time = (0:pnts - 1) / srate * 1000 - baseline;
% overlap = 0;
overlap = 100;

%% Classical Measures
CMwindow = struct('length', winLength, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', srate', 'baseline', 0);
rawCMconfig = struct('measures', [], 'time', time,  'freqRange', [],...
    'window', CMwindow, 'statistics', 0, 'nSurrogates', 100,...
    'fs', srate);
rawCMconfig.measures = {'COH', 'iCOH'};

cmIndexes = H_compute_CM_commandline(EEGdata, rawCMconfig);
cData = cell2mat(cmIndexes.COH.data);
icData = cell2mat(cmIndexes.iCOH.data);
cmIndexDimen = cmIndexes.COH.dimensions(6);
freqBands = {[1 4], [4 8], [8 12], [12 30], [30 46]};

clear cmIndexes rawCMconfig CMwindow

timeLength = size(cData,4);
for l = 1:bandsNum
    FrequencyBand = cell2mat(freqBands(l));
    disp(['CM processing for ', mat2str(FrequencyBand)]);
    f = cell2mat(cmIndexDimen);
    indexBand = find(floor(f)>= FrequencyBand(1,1)...
        &  floor(f)<= FrequencyBand(1,2));
    freqIndex = [indexBand(1) indexBand(end)];

    cohIndexes = zeros(chanNum, chanNum, 1, timeLength);
    icohIndexes = zeros(chanNum, chanNum, 1, timeLength);
    for k = 1:timeLength
        CoherenceMatrix = ones(chanNum,chanNum);
        iCoherenceMatrix = ones(chanNum,chanNum);
        for i = 1:chanNum
            for j=i+1:chanNum
                meas = reshape(...
                    cData(i,j,:,k),[size(cData(i,j,:,k),3) 1]);
                imeas = reshape(...
                    icData(i,j,:,k),[size(icData(i,j,:,k),3) 1]);
                tempVec = 0.5*log((1+meas)./(1-meas));
                itempVec = 0.5*log((1+imeas)./(1-imeas));
                CoherenceMatrix(i,j) = nanmean(...
                    tempVec(freqIndex(1):freqIndex(2)));
                CoherenceMatrix(i,j) = (exp(2*CoherenceMatrix(i,j))-1)...
                    ./(exp(2*CoherenceMatrix(i,j))+1);
                iCoherenceMatrix(i,j) = nanmean(...
                    itempVec(freqIndex(1):freqIndex(2)));
                iCoherenceMatrix(i,j) = (exp(2*iCoherenceMatrix(i,j))-1)...
                    ./(exp(2*iCoherenceMatrix(i,j))+1);
            end
        end
        cohIndexes(:,:,:,k) = (CoherenceMatrix + CoherenceMatrix') - 1;
        icohIndexes(:,:,:,k) = (iCoherenceMatrix + iCoherenceMatrix') - 1;
    end
    fBands.(cell2mat(fBandsCell(l))).COH.Data = mat2cell(...
        cohIndexes, chanNum, chanNum, 1, timeLength);
    fBands.(cell2mat(fBandsCell(l))).iCOH.Data = mat2cell(...
        icohIndexes, chanNum, chanNum, 1, timeLength);
end

clear cData icData cmIndexDimen freqBands FrequencyBand f indexBand...
    freqIndex cohIndexes icohIndexes CoherenceMatrix iCoherenceMatrix...
    meas imeas tempVec itempVec

%% Phase Synchronization
PSwindow = struct('length', winLength, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', srate', 'baseline', 0);

rawPSconfig = struct('measures', [], 'bandcenter', [],...
    'bandwidth', 4, 'fs', srate, 'method', 'ema', 'time', time,...
    'window', PSwindow, 'statistics', 0, 'nSurrogates', 100);
rawPSconfig.measures = {'PLV', 'PLI', 'RHO'};

bandcenterM = [2.5, 6, 10, 21, 38];
bandwidthM = [3, 4, 4, 18, 16];



disp('Running PS processing...');
for i = 1:bandsNum
    rawPSconfig.bandcenter = bandcenterM(i);
    rawPSconfig.bandwidth = bandwidthM(i);
    psIndexes = H_compute_PS_commandline(EEGdata, rawPSconfig);
    bandsProcessed = cell2mat(psIndexes.PLV.dimensions(6));
    disp(['Finished for ', mat2str(bandsProcessed)]);
    fBands.(cell2mat(fBandsCell(i))).PLI.Data = psIndexes.PLI.data;
    fBands.(cell2mat(fBandsCell(i))).PLV.Data = psIndexes.PLV.data;
    fBands.(cell2mat(fBandsCell(i))).RHO.Data = psIndexes.RHO.data;
end

clear psIndexes PSwindow rawPSconfig bandcenterM bandwidthM...
    bandsProcessed pliData plvData rhoData

%% Brain Connectivity Measures and Features Vector Creation

features = zeros(6,150);
featuresIndex = cell(6,150);

disp('Running BC measures processing...');
h = 0;
for i = 1:bandsNum
    for j = 1:measuresNum
        currentDat = cell2mat(fBands.(cell2mat(fBandsCell(i)))...
            .(cell2mat(measuresCell(j))).Data);
        currentBCM = zeros(bcMeasuresNum, timeLength);
        for k = 1:timeLength
            currentMat = weight_conversion(...
                reshape(currentDat(:,:,:,k),[chanNum chanNum]),'autofix');
            [~, currentBCM(1,k)] = modularity_und(currentMat);
            currentBCM(2,k) = transitivity_wu(currentMat);
            [currentBCM(3,k), currentBCM(4,k),...
                ~, currentBCM(5,k), currentBCM(6,k)] = charpath(...
                distance_wei(weight_conversion(currentMat,'lengths')));
        end
        for k = 2:bcMeasuresNum+1
            fBands.(cell2mat(fBandsCell(i))).(cell2mat(measuresCell(j)))...
                .(cell2mat(bcMeasuresCell(k))) = ...
                mat2cell(currentBCM(k-1,:), 1, timeLength);
            staMeaTemp = [mean(currentBCM(k-1,:), 2);...
                median(currentBCM(k-1,:), 2);...
                std(currentBCM(k-1,:), 0, 2);...
                iqr(currentBCM(k-1,:), 2);...
                skewness(currentBCM(k-1,:), 1, 2);...
                kurtosis(currentBCM(k-1,:),1,2)];
            features(:,k-1+h) = staMeaTemp;
            for ll = 1:staMeasuresNum
%                 fFeatures.(cell2mat(fBandsCell(i)))...
%                 .(cell2mat(measuresCell(j)))...
%                 .(cell2mat(bcMeasuresCell(k)))...
%                 .(cell2mat(staMeasuresCell(ll))) = staMeaTemp(ll);
                fIndexStr = strjoin({cell2mat(fBandsCell(i)),...
                    cell2mat(measuresCell(j)),...
                    cell2mat(bcMeasuresCell(k)),...
                    cell2mat(staMeasuresCell(ll))});
                featuresIndex(ll,k-1+h) = mat2cell(fIndexStr,...
                    size(fIndexStr, 1), size(fIndexStr, 2));
            end
        end
        h = h + 6;
    end
end

features = reshape(features, [1 900]);
featuresIndex = reshape(featuresIndex, [1 900]);

disp('Processing finished successfully');

end