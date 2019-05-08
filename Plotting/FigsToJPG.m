function FigsToJPG(directory)

if isfolder(directory)
    cd(directory)
    files = dir('*.fig');
    for i = 1:length(files)
        fig = openfig(files(i).name);
        name = erase(files(i).name,'.fig');
        export_fig(
        %saveas(fig,[name '.jpg']);
        close(fig)
    end
elseif isfile(directory)
    fig = openfig(directory);
    name = erase(directory,'.fig');
    saveas(fig,[name '.jpg']);
end