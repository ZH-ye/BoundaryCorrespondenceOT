% demo

load('data.mat')
addpath('./src');
addpath('./3rd');
addpath('./model');
addpath('./utils/')
%% find four corners for a polygon %%%%

name = 'cattle';
data = mat(:,:,1);
N = size(data,2); % number of points

cost_type = [0.5,1/16];
lambda = [0.05,0.05];

% below is a step by step demostration


%% initialization  %%%%%%%%%%%%%%%%%%
use_local_search = false;
[corners_init,~]=simp_sink(data,cost_type,lambda,use_local_search);

%% local search %%%%%%%%%%%%%%%%%%%%%
[~,E,chord_para]=local_search_mex(data,corners_init,cost_type,lambda,1:N); 

%% post process %%%%%%%%%%%%%%%%%%%%%

[corners,~] = project_to(data,chord_para,1:N);
% removed invalid ones
valid = all(diff(corners,1)>0,1);
corners = corners(:,valid);
E = E(valid);
% unique column
[corners,ia,~] = unique(corners','stable','rows');
corners = corners';
E = E(ia);
corners = sort(corners,1);

% all in one : [corners,~]=simp_sink(data,cost_type,lambda);

%% the final result %%%%%%%%%%%%%%%%%%%
c = corners(:,1);
cornerpoints = data(:,c);

%% plot %%%%

coeff = spline_parametrization(data,c);
spline_plot(name,data,c,coeff,'polygon')
spline_plot([name,' mu'],data,c,coeff,'uv','mu')
spline_plot([name,' log Jac'],data,c,coeff,'heatmap','log_jacobian')


