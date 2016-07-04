clear all

if ispc
    EEG = pop_biosig('H:\eeg_hei_soton\2.1.1min');
    EEG=pop_chanedit(EEG, 'lookup',['H:\\thesis\\eeglab13_5_4b\\plugins'...
        '\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp']);
else
    EEG = pop_biosig('/home/andres/repos/eeg_hei_soton/2.1.1min');
    EEG = pop_chanedit(EEG, 'lookup',['/home/andres/MATLAB/'...
        'eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/'...
        'standard-10-5-cap385.elp']);
end

%% Filtering and Epoching
EEG = pop_select( EEG,'nochannel',{'FPZ' 'AUX1' 'AUX2' 'AUX3' 'AUX4',...
    'AUX5' 'AUX6' 'AUX7' 'AUX8' 'PG1' 'PG2' 'A1' 'A2'});
notchlo = 49;
notchhi = 51;
bandlo = 0.5;
bandhi = 45;
epochlen = 5;
EEG = pop_eegfiltnew(EEG, notchlo, notchhi, 1690, 1, [], 0);
EEG = pop_eegfiltnew(EEG, bandlo, bandhi, 3380, 0, [], 0);
EEG.data(1,1:EEG.srate*epochlen:EEG.pnts) = 1;  % Set events
EEG = pop_chanevent(EEG, 1,'edge','leading','edgelen',0,'duration','on');
EEG = pop_epoch(EEG, {  }, [-1  2], 'newname', '1MIN data epochs',...
    'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-1000     0]);   % Check what Baseline removing is 
                                        % for, and proper Value!!!!
[EEG, rejectIndexes] = pop_eegthresh(EEG,1,1:19,-100,100,-1,1.998,0,0);
EEG = pop_rejepoch( EEG, rejectIndexes, 0);

% onemat = 0;
% if (onemat)
%     name = 'epochs.mat';
%     oData = double(EEG.data);
%     save(name,'oData');
% else    
%     for i = 1:EEG.trials
%         name = ['EEGLAB output/epochs/e', num2str(i),'.mat'];
%         oData = double(EEG.data(:,:,i));
%         save(name,'oData');
%     end
% end
% eeglab redraw

measuresCell = {'COH'; 'iCOH'; 'PLV'; 'PLI'; 'RHO'};
measuresDataCell = {[]; []; []; []; []};
measures = cell2struct(measuresDataCell, measuresCell, 1);

fBandsCell = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
for i = 1:size(fBandsCell,2)
    fBands.(cell2mat(fBandsCell(i))) = measures;
end

clear measures measuresDataCell measuresCell

%% Window Properties
% length = EEG.pnts*1000/EEG.srate; % Lenght in ms
length = 1000;
baseline = 0;
time = (0:EEG.pnts - 1) / EEG.srate * 1000 - baseline;
% overlap = 0;
overlap = 100;

%% Classical Measures
CMwindow = struct('length', length, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', EEG.srate', 'baseline', 0);
rawCMconfig = struct('measures', [], 'time', time,  'freqRange', [],...
    'window', CMwindow, 'statistics', 0, 'nSurrogates', 100,...
    'fs', EEG.srate);
rawCMconfig.measures = {'COH', 'iCOH'};

cmIndexes = H_compute_CM_commandline(EEG.data, rawCMconfig);
cData = cell2mat(cmIndexes.COH.data);
icData = cell2mat(cmIndexes.iCOH.data);
cmIndexDimen = cmIndexes.COH.dimensions(6);
freqBands = {[1 4], [4 8], [8 12], [12 30], [30 45]};

clear cmIndexes rawCMconfig CMwindow

for l = 1:size(freqBands,2)
    FrequencyBand = cell2mat(freqBands(l));
    disp(['CM processing for ', mat2str(FrequencyBand)]);
    f = cell2mat(cmIndexDimen);
    indexBand = find(floor(f)>= FrequencyBand(1,1)...
        &  floor(f)<= FrequencyBand(1,2));
    freqIndex = [indexBand(1) indexBand(end)];

    cohIndexes = zeros(size(cData,1), size(cData,2), 1, size(cData,4));
    icohIndexes = zeros(size(icData,1), size(icData,2), 1, size(icData,4));
    for k = 1:size(cData,4)
        CoherenceMatrix = ones(size(cData,1),size(cData,1));
        iCoherenceMatrix = ones(size(icData,1),size(icData,1));
        for i = 1:size(cData,1)
            for j=i+1:size(cData,2)
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
    fBands.(cell2mat(fBandsCell(l))).COH = mat2cell(...
        cohIndexes, size(cData,1), size(cData,2), 1, size(cData,4));
    fBands.(cell2mat(fBandsCell(l))).iCOH = mat2cell(...
        icohIndexes, size(cData,1), size(cData,2), 1, size(cData,4));
end

clear cData icData cmIndexDimen freqBands FrequencyBand f indexBand...
    freqIndex cohIndexes icohIndexes CoherenceMatrix iCoherenceMatrix...
    meas imeas tempVec itempVec

%% Phase Synchronization
PSwindow = struct('length', length, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', EEG.srate', 'baseline', 0);

rawPSconfig = struct('measures', [], 'bandcenter', [],...
    'bandwidth', 4, 'fs', EEG.srate, 'method', 'ema', 'time', time,...
    'window', PSwindow, 'statistics', 0, 'nSurrogates', 100);
rawPSconfig.measures = {'PLV', 'PLI', 'RHO'};

bandcenterM = [2.5, 6, 10, 21, 38];
bandwidthM = [3, 4, 4, 18, 16];

for i = 1:size(bandcenterM, 2)
    rawPSconfig.bandcenter = bandcenterM(i);
    rawPSconfig.bandwidth = bandwidthM(i);
    psIndexes = H_compute_PS_commandline(EEG.data, rawPSconfig);
    bandsProcessed = cell2mat(psIndexes.PLV.dimensions(6));
    disp(['PS processing for ', mat2str(bandsProcessed)]);    
    fBands.(cell2mat(fBandsCell(i))).PLI = psIndexes.PLI.data;
    fBands.(cell2mat(fBandsCell(i))).PLV = psIndexes.PLV.data;
    fBands.(cell2mat(fBandsCell(i))).RHO = psIndexes.RHO.data;
end

clear psIndexes PSwindow rawPSconfig bandcenterM bandwidthM...
    bandsProcessed pliData plvData rhoData

disp('Processing finished successfully');