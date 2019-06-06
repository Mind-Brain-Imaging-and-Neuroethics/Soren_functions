function align_folders(folders,suffixes)

for i = 1:length(folders)
    cd(folders{i})
    files = dir(['*' suffixes{i}]);
    names = extractfield(files,'name');
    subs{i} = erase(names,suffixes{i});
end

commonsubs = intersect(subs{1},subs{2});

if length(folders) > 2
    for i = 3:length(folders)
        commonsubs = intersect(commonsubs,subs{i});
    end
end

for i = 1:length(folders)
   cd(folders{i})
   mkdir('Unused')
   files = dir(['*' suffixes{i}]);
   for c = 1:length(files)
      subid = erase(files(c).name,suffixes{i});
      if ~strcmpi(subid,commonsubs)
          system(['mv ' files(c).name ' Unused'])
      end
   end
end
