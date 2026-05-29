function save_solution(solution, summary, out_file)
%SAVE_SOLUTION Persist optimal trajectories and summary metadata.

    if nargin < 3 || isempty(out_file)
        repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        out_file = fullfile(repo_root, 'results', 'latest_solution.mat');
    end

    out_dir = fileparts(out_file);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    save(out_file, 'solution', 'summary');
    fprintf('Saved results to %s\n', out_file);
end
