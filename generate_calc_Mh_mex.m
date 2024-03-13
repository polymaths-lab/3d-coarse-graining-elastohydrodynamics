%Run this script to generate calc_Mh_mex function.

% GENERATE_CALC_MH_MEX   Generate MEX-function calc_Mh_mex from calc_Mh.
% 
% Script generated from project 'calc_Mh.prj' on 13-Mar-2024.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'calc_Mh'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg calc_Mh -args ARGS{1}

