function [EEGout] = EEG_preprocess(filepath,subID)



eeglab

%load the data file
EEG = pop_fileio([filepath]);
EEG.setname = ['Grouping' num2str(subID)];
EEG = eeg_checkset( EEG );


%select the relevant data

EEG = pop_select(EEG,'point',[550000 950000]);


% %downsample to 250 Hz
% EEG = pop_resample( EEG, 250);
% EEG.setname=['Grouping' num2str(subID)];
% EEG = eeg_checkset( EEG );

%look up channel locations
EEG=pop_chanedit(EEG, 'lookup','/Users/Soren/Documents/MATLAB/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );

% %find all the listening trials (if I fix my triggers this might be less of a pain)
% [boundariesS1] = GroupTrials(EEG);
% boundaries = boundaries - 1;
% EEG = pop_select( EEG,'point',[EEG.events(2).latency - EEG.srate*30, EEG.events(listenbounds).latency + EEG.srate*30]);
% EEG.setname = ['Grouping' num2str(subID)];
% EEG = eeg_checkset( EEG );

%high pass filter
EEG = pop_eegfiltnew(EEG, [], 0.1, 826, true, [], 1);
EEG.setname=['Grouping' num2str(subID)];
EEG = eeg_checkset( EEG );

%ASR
EEG = clean_rawdata(EEG, 5, [-1], 0.8, -1, 5, 0.5);
EEG = eeg_checkset( EEG );

%PREP pipeline (remove line noise, re-reference with robust average
%reference)
EEG = pop_prepPipeline(EEG);
EEG.setname = ['Grouping' num2str(subID)];
EEG = eeg_checkset( EEG );

%ICA
EEG = pop_runica(EEG, 'extended',1,'interupt','on');