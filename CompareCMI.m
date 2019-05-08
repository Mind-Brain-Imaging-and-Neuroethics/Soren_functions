function [outmeasures] = CompareCMI(measurehandles,varargin)

listSubs = {'A00053990','A00054207','A00054387','A00054930','A00055065','A00059083','A00062919','A00054039','A00057092','A00058775'};
testSample = [];


for i = 1:length(listSubs)
   try
    load(['/Users/Soren/Desktop/qEEG_analysis/CMI_database/' listSubs{i} '/EEG/preprocessed/mat_format/gp_' listSubs{i} '001.mat'],'EEG');
   catch
          load(['/Users/Soren/Desktop/qEEG_analysis/CMI_database/' listSubs{i} '/EEG/preprocessed/mat_format/gip_' listSubs{i} '001.mat'],'EEG');
   end
   
   EEG = pop_eegfiltnew(EEG, 18, 36, 368, 0, [], 0);
   
   Signal = EEG.data;
   Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function
   SignalInfo = nbt_Info; %this initializes an Info Object
   SignalInfo.converted_sample_frequency = EEG.srate;
   AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, 0.5, 50, 4);
   EEG.data = AmplitudeEnvelope';
   
   for c = 1:length(measurehandles)
       outmeasures(i,:,c) = measurehandles{c}(EEG);
   end

   clear EEG
end

% for c = 1:length(measurehandles)
% meanMeasures(c) = nan_mean(outmeasures(:,:,c),2);
% end

