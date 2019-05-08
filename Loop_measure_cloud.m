function [outmeasures] = Loop_measure_cloud(measurehandles,filepath,varargin)

%standardized loop to apply a measure to an EEG data set

%measurehandle should be a function handle for the measure you want to
%apply: make sure to create a new function file for this measure that can
%take a whole EEG struct as input and output the data in the format
%subjects x channels
%filepath is the path to the file where you want to save the data



eeglab

if CheckInput(varargin,'Dir')
    directory_name = EasyParse(varargin,'Dir');
else
    directory_name = uigetdir;
end

cd(directory_name);

%%%%Let's select the relevant files, I am assuming here continuous resting-state files already clean
files = dir(fullfile(directory_name, '*.set'));

fileindex = find(~[files.isdir]);


%%%%Loop through all the files
startSub = 1;

if EasyParse(varargin,'Continue','true')
    load(filepath)
    startSub = size(outmeasures,1)+1;
end

for i = startSub:length(fileindex)
    filename = files(fileindex(i)).name;
    [PATH, NAME, EXT] = fileparts(filename);
    
    NAME = [NAME, EXT];
    
    disp(' ')
    disp(['Now processing subject ' num2str(i)])
    
    EEG = pop_loadset( 'filename', filename, 'filepath', directory_name);
    
    if nargin > 2
        if CheckInput(varargin,'Subsample')
            sampleSize = EasyParse(varargin,'Subsample');
            if EasyParse(varargin,'StartPoint','FirstAfterRest')
                for c = 4:length(EEG.event)-1
                    if strcmpi(EEG.event(c).type,'S  1') || strcmpi(EEG.event(c).type,'S 1')
                        startPoint = EEG.event(c+1).latency;
                        break;
                    end
                end
                
            elseif CheckInput(varargin,'StartPoint')
                startPoint = EasyParse(varargin,'StartPoint');
            else
                startPoint = randi(length(EEG.data-sampleSize));
            end
            
            disp(' ')
            disp(['Processing data from data point ' num2str(startPoint) ' to data point ' num2str(startPoint+sampleSize) '...'])
            
            EEG.data = EEG.data(:,startPoint:(startPoint+sampleSize));
        end
        
        if CheckInput(varargin,'Filter')
           if EasyParse(varargin,'Filter','Alpha')
               disp(' ')
               disp('Filtering 8 - 13 Hz...')
               EEG = pop_eegfiltnew(EEG, 8,13,826,0,[],0);
           end
        end
        
        
        if EasyParse(varargin,'AmpEn','true')
            disp(' ')
            disp('Getting amplitude envelope...')
            Signal = EEG.data;
            
            Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function
            
            SignalInfo = nbt_Info; %this initializes an Info Object
            SignalInfo.converted_sample_frequency = 500; %Sets the frequency to 500Hz.
            
            
            AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, 0.5, 50, 4);
            
            EEG.data = AmplitudeEnvelope';
        end
            
        if EasyParse(varargin,'IRASA','true')
            
            disp(' ')
            disp('Performing IRASA...')
            EEG = amri_sig_fractal(EEG.data',EEG.srate);
            
        end
        
    end
    

    
    
    for c = 1:length(measurehandles)
        outmeasures(i,:,c) = measurehandles{c}(EEG);
    end
    
    %save after every subject so you can continue later
    save(filepath,'outmeasures');
    
end
