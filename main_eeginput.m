clear global

subjects = {'2.1.1', '2.4.2', '2.5.1', '2.5.2', '2.2.1', '2.2.2', '2.3.2'};

for i = 1:size(subjects,2)
    sub = cell2mat(subjects(i));
    dotlessSub = ['sub' regexprep(sub, {'\.'},{''})];
    features.(dotlessSub) = eeg_hie_feaX(...
        ['H:\repos\eeg_hie_soton\RAWEEGDATA\' sub ' 1min']);
end