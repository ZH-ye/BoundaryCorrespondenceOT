%% simpleDomains
load('simpleDomains.mat')
addpath('./src');
addpath('./3rd');
addpath('./model');
addpath('./utils/')


name = 'domain';
for ii = 1:4
%% find four corners for simple domains %%%%
data = mat(:,:,ii);
N = size(data,2); % number of points

cost_type = [1,1/16];
lambda = [0.05,0.05];

[corners,~]=simp_sink(data,cost_type,lambda);

c = corners(:,1);
cornerpoints = data(:,c);

%% plot %%%%

coeff = spline_parametrization(data,c,0);
spline_plot(sprintf('%s - %d',name,ii),data,c,coeff,'uv','blue and red')


end