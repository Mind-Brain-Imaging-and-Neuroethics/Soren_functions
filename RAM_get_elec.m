function [elec,channum] = RAM_get_elec(filename)

contacts = jsonread(filename);
subid = fieldnames(contacts);
subid = subid{1};

fields = fieldnames(contacts.(subid).contacts);

for c = 1:length(fields)
    channum(c) = contacts.(subid).contacts.(fields{c}).channel;
    if isfield(contacts.(subid).contacts.(fields{c}).atlases,'mni')
        pos(c,1) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.mni.x);
        pos(c,2) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.mni.y);
        pos(c,3) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.mni.z);
    elseif isfield(contacts.(subid).contacts.(fields{c}).atlases,'tal')
        pos(c,1) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.tal.x);
        pos(c,2) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.tal.y);
        pos(c,3) = test_str2num(contacts.(subid).contacts.(fields{c}).atlases.tal.z);
    end
end
if ~isfield(contacts.(subid).contacts.(fields{c}).atlases,'mni') && isfield(contacts.(subid).contacts.(fields{c}).atlases,'tal')
    pos = tal2icbm_spm(pos);
end

[~,sorted] = sort(channum);

elec = struct;

elec.unit = 'mm';
elec.coordsys = 'mni';
elec.label = fields(sorted);
elec.elecpos = pos(sorted,:);
elec.chanpos = pos(sorted,:);
elec.tra = eye(length(elec.label));

end

function [output] = test_str2num(input)

if ischar(input) || isstr(input)
    output = str2num(input);
else
    output = input;
end

end
