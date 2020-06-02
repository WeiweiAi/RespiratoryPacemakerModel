warning off;
path_var=pwd;
%% Set up for the staliro
staliroDir='/nesi/nobackup/uoa00596/STALIRO/trunk'; 
cd(staliroDir)
setup_staliro;
%% Prepare the data need for the model
addpath(genpath(getenv('root_path'))); %%This will add $root_path and all subdirectories
cd(path_var)
dataDir=getenv('data_path'); % To save the simulation data
curdir=pwd;
path_var=getenv('working_path');
sid=str2num(getenv('SLURM_ARRAY_TASK_ID'))
model_name=getenv('model_name');
model = sprintf('%s%d',model_name,sid); % the model name
mdl=[path_var filesep model]; % the model with the full path.
varin=[7.5 0.1; 7.5 0.2; 7.5 0.3; 7.5 0.4;
	   7.0 0.1; 7.0 0.2; 7.0 0.3; 7.0 0.4;
	   6.5 0.1; 6.5 0.2; 6.5 0.3; 6.5 0.4;
	   6.0 0.1; 6.0 0.2; 6.0 0.3; 6.0 0.4;];
% Setup tempdir and cd into it
tmpDir = tempname;
mkdir(tmpDir);
cd(tmpDir);
path_var=pwd;

n_tests = 1000;
time = 2000;
init_cond = [];
Kp=varin(sid,2)
LRI=varin(sid,1)
input_range = [Kp 0.6; % Kp
               5 LRI]    % LRI             
cp_array=[1 1];

preds(1).str = 'lb'; 
preds(1).A = [-1];
preds(1).b = -95; %a constraint of the form Ax<=b; lb>=95
preds(1).loc = []; %If the control location vector is empty, then the predicate should hold in any location
% 
% preds(2).str = 'ub'; 
% preds(2).A = [1];
% preds(2).b = 100; % ub<=100
% preds(2).loc = []; %If the control location vector is empty, then the predicate should hold in any location

phi = '[](!lb-><>_[0,60] lb)'; 

disp(' ')
disp('Create an staliro_options object with the default options:')
opt = staliro_options();
opt.runs = 1;
opt.interpolationtype = {'const'};
opt.optimization_solver = 'SA_Taliro';
opt.falsification = 0; % stop at violation
opt.optim_params.n_tests = n_tests; % Maximum test number
opt
disp(' ')
disp('Running S-TaLiRo ...')
[results, history, opt] = staliro(model,init_cond,input_range,cp_array,phi,preds,time,opt);
filename=sprintf('simresults_req%d.mat',sid);
save(filename,'history','opt');
movefile([path_var filesep filename],dataDir);
% Switch all back to their original folder.
    clear functions
    close_system(model, 0);
    cd(curdir);
    rmdir(tmpDir,'s');
	
	