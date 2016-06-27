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

length = EEG.pnts*1000/EEG.srate; % Lenght in ms
baseline = 0;
time = (0:EEG.pnts - 1) / EEG.srate * 1000 - baseline;
overlap = 0;

CMwindow = struct('length', length, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', EEG.srate', 'baseline', 0);
rawCMconfig = struct('measures', [], 'time', time,  'freqRange', [],...
    'window', CMwindow, 'statistics', 0, 'nSurrogates', 100,...
    'fs', EEG.srate);
rawCMconfig.measures = {'COH', 'iCOH'};

cmIndexes = H_compute_CM_commandline(EEG.data, rawCMconfig);

PSwindow = struct('length', length, 'overlap', overlap,...
    'alignment', 'epoch', 'fs', EEG.srate', 'baseline', 0);

rawPSconfig = struct('measures', [], 'bandcenter', [10 20],...
    'bandwidth', 4, 'fs', EEG.srate, 'method', 'ema', 'time', time,...
    'window', PSwindow, 'statistics', 0, 'nSurrogates', 100);
rawPSconfig.measures = {'PLV', 'PLI', 'RHO'};

psIndexes = H_compute_PS_commandline(EEG.data, rawPSconfig);