function [outmeasures] = Loop_measure(measurehandles,filepath,varargin)

%Loop_measure: a standardized loop to apply a measure to an EEG data set

%Required input parameters:
%
%measurehandle should be a function handle for the measure you want to
%apply: make sure to create a new function file for this measure that can
%take a whole EEG struct as input and output the data in the format
%subjects x channels
%
%filepath is the path to the file where you want to save the data
%
%Optional key-value pairs:
%'Format','.mat': changes the script so it runs on .mat files rather than
%EEG data sets - I used this for the LORETA stuff
%'Format','Fieldtrip': uses .mat files which contain a fieldtrip structure
%
%'Dir',string: you can give the script a directory in string form if you
%don't want to select it in the popup window
%
%'Continue','true': the script will load the file from the filepath and
%continue from the last subject that was saved
%
%'Filter','Alpha': filters the dataset from 8-13 Hz before continuing
%
%'Spectrum','true': calculates the power spectrum and add the power and
%frequencies to the EEG struct
%
%'AmpEn','true': calculates the amplitude envelope and changes EEG.data to
%the amplitude envelope. All input measures are then calculated on this
%amplitude envelope
%
%'IRASA','true': performs the IRASA method, and then changes the EEG struct
%to the IRASA spect struct - this means you can only use measurehandle
%inputs that take this IRASA spectrum as input
%
%'Subsample',integer: takes a subsample of the data to calculate measures
%on. Useful for running horizontal analysis on task data. Uses n samples
%
%'StartPoint',integer: sets the start point in the data from which the n
%samples in subsample are taken

%Outputs
%
%outmeasures: 3-D array of measures.
%   row: the subject (or file index where the subject was found)
%   column: the channel
%   3rd dimension: the measures being applied (in the order you input them)
%
%if surrogating, outmeasures will instead be a cell array of 3-D arrays
%with the following structure:
%   cell: the subject
%   row: the channel
%   column: the value for surrogate n
%   3rd dimension: the measures
%
%
%Make sure you add the folder 'FunctionWrappers' to your Matlab path. Any
%function in there can be used with this function



if ~CheckInput(varargin,'Format')
    eeglab
end

if CheckInput(varargin,'Dir')
    directory_name = EasyParse(varargin,'Dir');
else
    directory_name = uigetdir;
end

cd(directory_name);

if EasyParse(varargin,'Format','.mat')
    files = dir(fullfile(directory_name, '*.mat'));
elseif EasyParse(varargin,'Format','Fieldtrip')
    if CheckInput(varargin,'VarString')
            files = dir(fullfile(directory_name,EasyParse(varargin,'VarString')));
    else
    files = dir(fullfile(directory_name, '*.mat'));
    end
else
    files = dir(fullfile(directory_name, '*.set'));
end

fileindex = find(~[files.isdir]);


%%%%Loop through all the files
startSub = 1;

if EasyParse(varargin,'Continue','true')
    load(filepath)
    if ~CheckInput(varargin,'Surrogate')
        startSub = size(outmeasures,1)+1;
    else
        startSub = length(outmeasures)+1;
    end
end

argsin = varargin;

for i = startSub:length(fileindex)
    filename = files(fileindex(i)).name;
    
    disp(' ')
    disp(['Now processing subject ' num2str(i)])
    
    if EasyParse(varargin,'Format','.mat')
        EEG = struct;
        load(filename);
        EEG.data = thingy;
        EEG.srate = 500; %set this parameter manually each time
        EEG.nbchan = 5; %set this parameter manually each time
        EEG.trials = 1;
    elseif EasyParse(varargin,'Format','Fieldtrip')
        EEG = struct;
        allvars = load(filename);
        names = fieldnames(allvars);
        for c = 1:length(names)
           if isstruct(allvars.(names{c}))
              data = allvars.(names{c});
              clear allvars
              break
           end
        end
        clear names
        EEG.data = cat(2,data.trial{:});
        EEG.srate = data.fsample;
        EEG.nbchan = size(EEG.data,1);
        EEG.trials = 1;
        clear data
    else
        EEG = pop_loadset( 'filename', filename, 'filepath', directory_name);
    end
    
    if nargin > 2
        
        if CheckInput(varargin,'Surrogate')
            disp(['Analyzing with ' num2str(EasyParse(varargin,'Surrogate')) ' surrogates...'])
        end
        
        if CheckInput(varargin,'Subsample') && ~CheckInput(varargin,'MultiSample')
            EEG = SubSample(EEG,argsin)
        end
        
        if CheckInput(varargin,'Filter')
            if EasyParse(varargin,'Filter','Alpha')
                disp(' ')
                disp('Filtering 8 - 13 Hz...')
                EEG = pop_eegfiltnew(EEG, 8,13,826,0,[],0);
            elseif EasyParse(varargin,'Filter','HighBeta')
                disp(' ')
                disp('Filtering 18 - 36 Hz...')
                EEG = pop_eegfiltnew(EEG, 18, 36, 368, 0, [], 0); %only works for data sampled at 500 Hz
            end
        end
        
        if EasyParse(varargin,'Spectrum','true')
            disp(' ')
            disp('Calculating power spectrum...')
            
            [spectra,freqs] = spectopo(EEG.data,0,EEG.srate,'plot','off');
            
            EEG.spectra = spectra;
            EEG.freqs = freqs;
            
        end
        
        if CheckInput(varargin,'AmpEn')
            disp(' ')
            disp('Getting amplitude envelope...')
            Signal = EEG.data;
            
            Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function
            
            SignalInfo = nbt_Info; %this initializes an Info Object
            SignalInfo.converted_sample_frequency = EEG.srate;
            
            bpfreqs = EasyParse(varargin,'AmpEn');
            AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, bpfreqs(1), bpfreqs(2), 2*(1/bpfreqs(1)));
            
            EEG.data = AmplitudeEnvelope';
        end
        
        if EasyParse(varargin,'IRASA','true')
            
            disp(' ')
            disp('Performing IRASA...')
            EEG = amri_sig_fractal(EEG.data',EEG.srate);
            
        end
        
    end
    
    
    
    if ~CheckInput(varargin,'Surrogate') && ~CheckInput(varargin,'Multisample')
        for c = 1:length(measurehandles)
            outmeasures(i,:,c) = measurehandles{c}(EEG);
        end
    elseif CheckInput(varargin,'Multisample')
        disp(['Taking ' num2str(EasyParse(varargin,'Multisample')) ' samples...'])
        for q = 1:EasyParse(varargin,'Multisample')
            newEEG = SubSample(EEG,varargin);
            fprintf([num2str(q) ' '])
            for c = 1:length(measurehandles)
                outmeasures(i,:,q,c) = measurehandles{c}(newEEG);
            end
        end
    else
        for c = 1:EEG.nbchan
            if EasyParse(varargin,'SurrogateMethod','AAFT')
                disp(' ')
                disp(['Creating surrogates for channel ' num2str(c)])
                tmp = AAFT(EEG.data(c,:),EasyParse(varargin,'Surrogate'));
                EEG.data = tmp';
                EEG.nbchan = EasyParse(varargin,'Surrogate');
                disp(' ')
            else
                disp(' ')
                disp(['Creating surrogates for channel ' num2str(c)])
                tmp = zeros(EasyParse(varargin,'Surrogate'),length(EEG.data));
                for q = 1:EasyParse(varargin,'Surrogate')
                    tmp(q,:) = EEG.data(c,randperm(length(EEG.data)));
                end
                EEG.data = tmp;
                EEG.nbchan = EasyParse(varargin,'Surrogate');
            end
            for cc = 1:length(measurehandles)
                outmeasures{i}(c,:,cc) = measurehandles{cc}(EEG);
            end
        end
    end
    
    %save after every subject so you can continue later
    try
        save(filepath,'outmeasures','files');
    catch
        warning('Saving failed')
    end
    
end
end


function [EEG] = SubSample(EEG,argsin)
sampleSize = EasyParse(argsin,'Subsample');

if CheckInput(argsin,'StartBlock') %startblock is specifically for Ivar's data - specify a block to start from
    
    startIndex = [];
    for j = 1:length(EEG.event)-1
        if strcmpi(EEG.event(j).type,'S  1') || strcmpi(EEG.event(j).type,'S 1')
            startIndex = [startIndex EEG.event(j+1).latency];
        end
    end
    
    startPoint = startIndex(EasyParse(argsin,'StartBlock'));
    
elseif CheckInput(argsin,'StartPoint') %specify a specific latency to start at
    startPoint = EasyParse(argsin,'StartPoint');
else %if no startpoint specified, pick a random one
    startPoint = randi(length(EEG.data)-sampleSize);
end

disp(' ')
disp(['Processing data from data point ' num2str(startPoint) ' to data point ' num2str(startPoint+sampleSize) '...'])

if startPoint+sampleSize < length(EEG.data)
    EEG.data = EEG.data(:,startPoint:(startPoint+sampleSize));
else
    warning('Not enough samples to meet sample size requirement')
    disp(' ')
    disp(['Continuing with ' num2str(length(EEG.data)-startPoint) ' samples...'])
    disp([num2str(startPoint+sampleSize-length(EEG.data)) ' samples missing...'])
    EEG.data = EEG.data(:,startPoint:end);
end
end