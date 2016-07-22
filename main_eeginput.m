clear all
close all

subjects = {'2.1.1', '2.4.2', '2.5.1', '2.5.2', '2.2.1', '2.2.2', '2.3.2'};

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

for i = 1:size(subjects,2)
    sub = cell2mat(subjects(i));
    dotlessSub = ['sub' regexprep(sub, {'\.'},{''})];
    [features.(dotlessSub).data, features.(dotlessSub).index] = eeg_hie_feaX(...
        ['H:\repos\eeg_hie_soton\RAWEEGDATA\' sub ' 1min']);
    display('All subjects were successfully analysed.');
end