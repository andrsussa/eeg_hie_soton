clear
if ispc
    EEG = pop_biosig('H:\eeg_hei_soton\2.1.1min');
    EEG=pop_chanedit(EEG, 'lookup','H:\\thesis\\eeglab13_5_4b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
else
    EEG = pop_biosig('/home/andres/repos/eeg_hei_soton/2.1.1min');
    EEG = pop_chanedit(EEG, 'lookup','/home/andres/MATLAB/eeglab13_5_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
end
EEG = pop_select( EEG,'nochannel',{'FPZ' 'AUX1' 'AUX2' 'AUX3' 'AUX4' 'AUX5' 'AUX6' 'AUX7' 'AUX8' 'PG1' 'PG2' 'A1' 'A2'});
notchlo = 49;
notchhi = 51;
bandlo = 0.5;
bandhi = 45;
epochlen = 5;
EEG = pop_eegfiltnew(EEG, notchlo, notchhi, 1690, 1, [], 0);
EEG = pop_eegfiltnew(EEG, bandlo, bandhi, 3380, 0, [], 0);
EEG.data(1,1:EEG.srate*epochlen:EEG.pnts) = 1;
EEG = pop_chanevent(EEG, 1,'edge','leading','edgelen',0,'duration','on');
EEG = pop_epoch(EEG, {  }, [-1  2], 'newname', '1MIN data epochs', 'epochinfo', 'yes');
EEG = pop_rmbase(EEG, [-1000     0]); % Check what Baseline removing is for and proper Value!!!!
[EEG, Indexes] = pop_eegthresh(EEG,1,1:19,-100,100,-1,1.998,0,0);
EEG = pop_rejepoch( EEG, Indexes, 0);

onemat = 1;
if (onemat)
    name = 'epochs.mat';
    var = double(EEG.data);
    save(name,'var');
else    
    for i = 1:EEG.trials
        name = ['epochs/e', num2str(i),'.mat'];
        var = double(EEG.data(:,:,i));
        save(name,'var');
    end
end
% eeglab redraw