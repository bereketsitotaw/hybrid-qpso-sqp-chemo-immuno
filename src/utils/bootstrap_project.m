function repo_root = bootstrap_project()
%BOOTSTRAP_PROJECT Add repo paths and return repository root.
%
% If called before paths are set up, locates the repo from this file's path.
% Entry scripts should prefer init_project() in scripts/.

    if exist('init_project', 'file')
        repo_root = init_project();
        return;
    end

    this_file = mfilename('fullpath');
    repo_root = fileparts(fileparts(fileparts(this_file)));
    addpath(fullfile(repo_root, 'src', 'utils'));
    addpath(fullfile(repo_root, 'src', 'chemo_immuno'));
    addpath(fullfile(repo_root, 'configs'));
    addpath(fullfile(repo_root, 'tests'));
    addpath(fullfile(repo_root, 'scripts'));
    setup_project_paths();
end
