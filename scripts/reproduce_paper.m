%REPRODUCE_PAPER One-command reproduction of Section 5 case studies + Table 3.
%
% Usage (from repo root in MATLAB):
%   run('scripts/reproduce_paper.m')
%
% Programmatic options:
%   reproduce_paper_pipeline('freeze_reference', true)

clear; clc; close all;

addpath(fileparts(mfilename('fullpath')));
init_project();
reproduce_paper_pipeline();
