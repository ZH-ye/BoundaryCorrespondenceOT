%% more examples

load('MoreExamples.mat')
addpath('./src');
addpath('./3rd');
addpath('./model');
addpath('./utils/');

name = '';

data = mat(:,:,3);
N = size(data,2); % number of points

cost_type = [0.5,1/16];
lambda = [0.05,0.01];

[corners,~]=simp_sink(data,cost_type,lambda);

c = corners(:,1);
cornerpoints = data(:,c);

%% plot %%%%

coeff = spline_parametrization(data,c,10);
spline_plot(name,data,c,coeff,'uv','blue and red')