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
for c = 1:length(files)
    if contains(files(c).name,experiment)
    chan = str2num(char(extractAfter(files(c).name,'.')));
    if ismember(chan,channum)
        fid = fopen(fullfile(sessiondir,'noreref',files(c).name));
        cont_data.trial{1}(find(chan == sortchans),:) = fread(fid,sources.(sourcefield).data_format);
    end
    end
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




