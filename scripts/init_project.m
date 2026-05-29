function repo_root = init_project()
%INIT_PROJECT Add repository paths. Call at the top of every entry script.
%
% Resolves the repo root from this file's location (scripts/init_project.m),
% so it works regardless of MATLAB's current working directory.

    this_file = mfilename('fullpath');
    if isempty(this_file)
        error('init_project:NotFound', ...
            ['init_project.m must be on the MATLAB path. ', ...
             'From the repo root run: run(''scripts/run_case_a.m'')']);
    end

    script_dir = fileparts(this_file);
    repo_root = fileparts(script_dir);

    addpath(fullfile(repo_root, 'src', 'utils'));
    addpath(fullfile(repo_root, 'src', 'chemo_immuno'));
    addpath(fullfile(repo_root, 'configs'));
    addpath(fullfile(repo_root, 'tests'));
    addpath(fullfile(repo_root, 'scripts'));

    if exist('setup_project_paths', 'file')
        setup_project_paths();
    end
end
