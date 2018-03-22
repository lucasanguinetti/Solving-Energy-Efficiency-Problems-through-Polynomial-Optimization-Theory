%This Matlab script can be used to generate Figure 2 in the article:
%
%Andrea Pizzo, Alessio Zappone and Luca Sanguinetti, "Solving Energy Efficiency Problems
%through Polynomial Optimization Theory," IEEE Signal Processing Letters, Submitted to.
%
%This is version 1.0 (Last edited: 2018-22-03)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%

%Initialization
clear;
close all;
clc;
yalmip('clear')

%Order of the SOS relaxations (A tighter relaxation can be obtained by using a higher order relaxation)
SOSorderRelaxation_vec = 3:7;
maxiter = 10;

error_matrix = NaN*ones(length(SOSorderRelaxation_vec),maxiter);
niter_vec = 1:maxiter;

for ell=1:length(SOSorderRelaxation_vec)
    
    % run Algorithm 1
    SOSorderRelaxation = SOSorderRelaxation_vec(ell);
    yalmip('clear')
    [r_vec,xstar_vec,ystar_vec,lambda_vec,xopt_numerical,yopt_numerical,ropt_numerical] = Fig2simulation(SOSorderRelaxation, maxiter);
    
    % save objective values vs l vs # of iterations
    error_vec = (abs(ropt_numerical-r_vec*1e6))/ropt_numerical;
    error_matrix(ell,:) = error_vec;
    
end

%% plot
% create new figure
figure(1);

% Cell array of colors
Colarray = [0.8980, 0.1686, 0.3137, 0.5882, 0.2941, 0, 0, 0, .5, 1, 0.6, 0, 0.3, 1, 0.1, 0.2941, 0.33, 0.5];
% Cell array of markers
Markarray = ['.','*','o','d','p','x','>','<'];
for ell=1:length(SOSorderRelaxation_vec)

    semilogy(niter_vec,abs(smooth(error_matrix(ell,:),'sgolay',4)),'color',Colarray(3*(ell-1)+1:3*ell),'Marker',Markarray(ell),'MarkerSize',14,'Linewidth',3)
    legendInfo{ell} = ['data' num2str(ell)]; % or whatever is appropriate
    hold on;
end   
legend(legendInfo)
xticks(niter_vec)
xlabel('xaxis');
ylabel('yaxis');
set(gca,'Fontsize',20)
box on; grid on;
