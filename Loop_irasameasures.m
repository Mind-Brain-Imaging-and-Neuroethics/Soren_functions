function [outmeasures] = Loop_measure(measurehandles,filepath,AmpEn)

%standardized loop to apply a measure to an EEG data set

%measurehandle should be a function handle for the measure you want to
%apply: make sure to create a new function file for this measure that can
%take a whole EEG struct as input and output the data in the format
%subjects x channels
%filepath is the path to the file where you want to save the data



eeglab
directory_name = uigetdir;
cd(directory_name);

%%%%Let's select the relevant files, I am assuming here continuous resting-state files already clean
files = dir(fullfile(directory_name, '*.set'));

fileindex = find(~[files.isdir]);


%%%%Loop through all the files
for i = 1:length(fileindex)
    disp(i)
    filename = files(fileindex(i)).name;
    [PATH, NAME, EXT] = fileparts(filename);
    
    NAME = [NAME, EXT];
    
    disp(['Now processing subject ' num2str(i)])
    
    EEG = pop_loadset( 'filename', filename, 'filepath', directory_name);
    
    if nargin > 2 && AmpEn
        Signal = EEG.data;
        
        Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function
        
        SignalInfo = nbt_Info; %this initializes an Info Object
        SignalInfo.converted_sample_frequency = 500; %Sets the frequency to 500Hz.
        
        
        AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, 0.5, 50, 4);
        
        EEG.data = AmplitudeEnvelope';
    end
    
    
    for c = 1:length(measurehandles)
        outmeasures(i,:,c) = measurehandles{c}(EEG);
    end
    
    
    
end

save(filepath,'outmeasures');