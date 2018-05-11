for i = 1:size(dep,1)
    mkdir(fileparts(strrep(dep{i}, '/home/cnoslab-ra/Applications/spm12', pwd)))
    copyfile(dep{i}, strrep(dep{i}, '/home/cnoslab-ra/Applications/spm12', pwd))
end