function setup_project_paths()
%SETUP_PROJECT_PATHS Add src directories to the MATLAB path for this repo.

    this_file = mfilename('fullpath');
    repo_root = fileparts(fileparts(fileparts(this_file)));
    addpath(fullfile(repo_root, 'src', 'chemo_immuno'));
    addpath(fullfile(repo_root, 'src', 'utils'));
    addpath(fullfile(repo_root, 'configs'));
    addpath(fullfile(repo_root, 'tests'));
    addpath(fullfile(repo_root, 'scripts'));
end
