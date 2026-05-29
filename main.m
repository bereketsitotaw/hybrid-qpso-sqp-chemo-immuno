function main(case_id)
%MAIN Hybrid QPSO-SQP chemo-immunotherapy optimal control (Cases A/B/C).
%
% Usage (cd to repo root in MATLAB, then):
%   main          % run all three cases
%   main('a')     % Case A only (balanced objectives)
%   main('b')     % Case B only (asymmetric regularization)
%   main('c')     % Case C only (minimum-time limit)
%
% Paper: doi:10.3934/dcdsb.2026072

    if nargin < 1 || isempty(case_id)
        case_id = 'all';
    end

    clc;
    close all;
    init_project();
    run_cases(case_id);
end
