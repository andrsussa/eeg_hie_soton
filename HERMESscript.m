clear all
load('/home/andres/repos/eeg_hei_soton/EEGLAB output/cTEST.mat')
fs = 512;
window = struct('length', 3000, 'overlap', 0, 'alignment', 'epoch', 'baseline', 0, 'fs', fs);
rawconfig = struct('window', window, 'fs', fs, 'statistics', 0, 'nSurrogates', 100, 'time', 0:1000/fs:3000, 'freqRange', []);
rawconfig.measures = {'COH', 'iCOH'};
TESTindexes = H_compute_CM_commandline(cTEST(:,:,1), rawconfig);