function main(case_id)
%MAIN Hybrid QPSO-SQP chemo-immunotherapy optimal control (Cases A/B/C).
%
% Usage (cd to repo root in MATLAB, then):
%   main              % run all three paper cases
%   main('a')         % Case A only (balanced objectives)
%   main('b')         % Case B only (asymmetric regularization)
%   main('c')         % Case C only (minimum-time limit)
%   main('test')      % MATLAB smoke tests (no full paper solve)
%   main('smoke')     % tiny non-paper plumbing solve (Case A)
%   main('validate')  % Table 3 pin self-check (no solve)
%
% Paper: doi:10.3934/dcdsb.2026072

    if nargin < 1 || isempty(case_id)
        case_id = 'all';
    end

    clc;
    close all;
    init_project();

    key = lower(char(case_id));
    switch key
        case {'test', 'tests'}
            here = fileparts(mfilename('fullpath'));
            run(fullfile(here, 'tests', 'run_smoke_tests.m'));
        case 'smoke'
            run_cases('smoke');
        case {'validate', 'check'}
            report = validate_reference_pin('verbose', true);
            if ~report.all_passed
                error('main:ValidateFailed', 'Reference pin self-check failed.');
            end
        otherwise
            run_cases(case_id);
    end
end
