%% Run both TCE bioremediation scenarios
% Run this file from the repository root. Generated outputs are written to
% results/generated/method1 and results/generated/method2.

clear;
close all;
clc;

repo_root = fileparts(mfilename('fullpath'));
addpath(fullfile(repo_root, 'src'));

output_root = fullfile(repo_root, 'results', 'generated');
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

fprintf('Running Method 1: distributed biostimulation/bioaugmentation...\n');
method1 = run_method1_distributed(output_root); %#ok<NASGU>

fprintf('Running Method 2: permeable reactive biobarrier...\n');
method2 = run_method2_prb(output_root); %#ok<NASGU>

fprintf('Finished. Outputs are available in:\n%s\n', output_root);
