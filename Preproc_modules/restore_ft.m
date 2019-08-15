function ft = restore_ft(origft,ft)

fields = {'trialinfo','sampleinfo','grad','elec'};

for c = 1:length(fields)
   if isfield(origft,fields{c})
      ft.(fields{c}) = origft.(fields{c}); 
   end
end