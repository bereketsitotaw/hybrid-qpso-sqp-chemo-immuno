function repo_root = init_project()
%INIT_PROJECT Add repo folders to the MATLAB path.

    repo_root = fileparts(mfilename('fullpath'));
    addpath(fullfile(repo_root, 'src', 'chemo_immuno'));
    addpath(fullfile(repo_root, 'src', 'utils'));
    addpath(fullfile(repo_root, 'configs'));
end
