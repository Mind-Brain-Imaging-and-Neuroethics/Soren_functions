function [cont_data,latencies] = RAM_loaddata(subid,experiment,session,eventcode)

dirname = pwd;

[elec,channum] = RAM_get_elec(fullfile(dirname,subid,'localizations','0','montages','0','neuroradiology','current_processed','contacts.json'));

sessiondir = fullfile(dirname,subid,'experiments',experiment,'sessions',session,'ephys','current_processed');

files = dir(fullfile(sessiondir,'noreref','*'));
files = files(find(~cell2mat(extractfield(files,'isdir'))));
sources = jsonread(fullfile(sessiondir,'sources.json'));
sourcefield = fieldnames(sources);
sourcefield = sourcefield{1};

sortchans = sort(channum);
rmindex = [];
for c = 1:length(files)
    fprintf([num2str(c) ' '])
    if contains(files(c).name,experiment)
        filechans(c) = str2num(char(extractAfter(files(c).name,'.')));
        if ~ismember(filechans(c),filechans(1:c-1)) % avoid duplicates
            if ismember(filechans(c),channum) % make sure you have corresponding electrode information
                fid = fopen(fullfile(sessiondir,'noreref',files(c).name));
                cont_data.trial{1}(find(filechans(c) == sortchans),:) = fread(fid,sources.(sourcefield).data_format);
            else
                filechans(c) = NaN;
            end
        else
            filechans(c) = NaN;
        end
    end
end

filechans(find(isnan(filechans))) = [];

if size(cont_data.trial{1},1) ~= length(elec.label) %delete electrodes with no corresponding data
    rmindex = find(~ismember(sortchans,filechans));
    elec.label(rmindex) = [];
    elec.elecpos(rmindex,:) = [];
    elec.chanpos(rmindex,:) = [];
    elec.tra(rmindex,:) = [];
    elec.tra(:,rmindex) = [];
end

srate = sources.(sourcefield).sample_rate;
cont_data.time = {linspace(0,length(cont_data.trial{1})/srate,length(cont_data.trial{1}))};
cont_data.elec = elec;
cont_data.label = elec.label;
cont_data.fsample = srate;

events = jsonread(fullfile(subid,'experiments',experiment,'sessions',session,'behavioral','current_processed','task_events.json'));

eventtypes = extractfield(events,'type');
eventindx = find(strcmpi(eventcode,eventtypes));
latencies = extractfield(events,'eegoffset');
latencies = latencies(eventindx);




