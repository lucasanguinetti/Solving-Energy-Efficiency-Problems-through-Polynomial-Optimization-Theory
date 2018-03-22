function [ r_poly_vec, xstar_poly_vec, ystar_poly_vec, lambda_poly_vec, xopt_numerical, yopt_numerical, ropt_numerical ] = Fig2simulation( SOSorderRelaxation, maxiter )

%Retrieve the optimizers $x^star = [K^star, M^star]$ of Problem (1) by
%using both the iterative Algorithm 1 (red square) and exhaustive search (black triangle).
%At each iteration of Algorithm 1, the optimal is computed by solving Problem (21) 
%provided in Section II of the article: 
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
%% 
%INPUT: 
%SOSorderRelaxation     = order of the SOS relaxation (indicated as l in the article cited above)
%maxiter                = maximum number of iteration allowed for Algorithm 1 
% 
%OUTPUT: 
%r_poly_vec             = Global Optimum objective computed using Algorithm 1 
%xstar_poly_vec         = Golbal Optimum point (K^\star) computed using Algorithm 1 
%ystar_poly_vec         = Golbal Optimum point (M^\star) computed using Algorithm 1 
%lambda_poly_vec        = Dinkelbach's parameter of Algorithm 1 
%xopt_numerical         = Golbal Optimum point (K^\star) computed using Exhaustive search
%yopt_numerical         = Golbal Optimum point (M^\star) computed using Exhaustive search
%ropt_numerical         = Global Optimum objective computed using Exhaustive search


%% Load simulation parameters
Fig2define;

% *********************************************
%% Compute the optimal solution to the EE maximization problem of [8] by using Exhaustive search
%
%[8] E. Björnson, L. Sanguinetti, and M. Kountouris, “Deploying dense
%networks for maximal energy efficiency: Small cells meet massive
%MIMO,” IEEE Journal on Selected Areas in Communications, vol. 34,
%no. 4, pp. 832–847, April 2016.
%
% *********************************************

%Placeholders for storing of simulation results
EEmatrix = zeros(Mmax,Kmax); %EE for different lambda and gamma values (using theoretical formulas)
zetaOptimal = zeros(Mmax,Kmax); %Store optimal beta for each point using the theoretical formulas
ASEmatrix = zeros(Mmax,Kmax);
AECmatrix = zeros(Mmax,Kmax);
kvec = [1:Kmax];
mvec = [1:Mmax];

%Go through all number of users
for kindex = 1:Kmax
    
    k = kvec(kindex);
    %Go through all number of antennas
    for mindex = 1:Mmax
        
        m = mvec(mindex);
        %Compute B1bar and B2bar from Eq. (18) and Eq. (19), respectively
        B1 = k*(4/(alphaval-2)^2 + 1/(alphaval-1) + 2/(alphaval-2)) + m/(alphaval-1) + 2/(alphaval-2)/SNR;
        B2 = k*(1 + 2/(alphaval-2))*(1 + 1/SNR) + (1 + 1/SNR)/SNR;
        
        %Check if the two constraints in Eq. (23) are satisfied
        if B1*gammaval/(m-B2*gammaval)>=1 && (k*B1*gammaval/(m-B2*gammaval))<=tau
            
            %Compute and store the optimal zeta
            zeta = B1*gammaval/(m-B2*gammaval); % optimal pilot reuse factor
            zetaOptimal(mindex,kindex) = zeta;
            
            %Compute the objective function in Eq. (23) and store the EE
            ASE = k *(1-(B1*gammaval / (m-B2*gammaval))*k/tau) * log2(1+ gammaval);
            AEC = C0cal + C1barcal*k + D0cal*m + D1cal*m*k + k *(1-(B1*gammaval / (m-B2*gammaval))*k/tau)*(log2(1+ gammaval)*Acal*Bw + ucal);
            %
            EEmatrix(mindex,kindex) = Bw*ASE/AEC/1e6;
            
            % save ASE and AEC
            ASEmatrix(mindex,kindex) = ASE;
            AECmatrix(mindex,kindex) = AEC;
            
        end
        
    end
    
end

kvec_scaled = kvec;
mvec_scaled = mvec;

% Compute the EE-optimal of [8] by using Exhaustive search
[EEmaxM,optM] = max(EEmatrix(:,:),[],1);
[EEmax,optK] = max(EEmaxM);

%Display optimizers
Mopt = mvec_scaled(optM(optK));
Kopt = kvec_scaled(optK);

% %Density of the lines that are used in the 3d plot to make it easier to
% %see the shape
% gridDensity = 10;
% 
% %Plot Figure 3 from the paper
% figure;
% hold on; box on; grid on;
% surface(kvec_scaled,mvec_scaled,EEmatrix,'EdgeColor','none');
% colormap(autumn);
% zlim([0 6])
% 
% view([-148 25]);
% 
% xlabel('Number of UEs (K)')
% ylabel('Number of BS antennas (M)');
% zlabel('Energy efficiency [Mbit/Joule]');
% title('Numerical optimization')
% set(gca,'Fontsize',20)
%
% %Plot lines on top of the 3d surface, to make it easier to see the shape
% for m = [1 gridDensity:gridDensity:Mmax]
%     plot3(kvec_scaled,mvec_scaled(m)*ones(1,Kmax),EEmatrix(m,:),'k-');
% end
% 
% for k = [1 gridDensity:gridDensity:Kmax]
%     plot3(kvec_scaled(k)*ones(1,Mmax),mvec_scaled,EEmatrix(:,k),'k-');
% end

%Plot the optimal solution to the EE maximization problem
% plot3(kvec_scaled(optK),mvec_scaled(optM(optK)),EEmax,'kv','MarkerSize',20,'MarkerFaceColor','black');
% legend(['\gamma = ' num2str(gammaval)])


% *********************************************
%% Compute the optimal solution to the equivalent FPP (1) by using Exhaustive search
% *********************************************
xrange = [1:Kmax];
yrange = [1:Mmax];
[xgrid, ygrid] = meshgrid(xrange,yrange);
obj_fun = zeros(size(xgrid));
for ii=1:length(xrange)
    x = xrange(ii);
    for jj=1:length(yrange)
        y = yrange(jj);
        if (h00+h10*x+h01*y)>=0 && (q00+q10*x+q01*y+q20*x^2+q11*x*y)>=0
            obj_fun(jj,ii) = (f10*x + f11*x*y + f20*x^2 + f21*x^2*y + f30*x^3)...
                /(g00 + g10*x + g01*y + g20*x^2 + g02*y^2 + g11*x*y + g21*x^2*y + g12*x*y^2 + g30*x^3);
        end
    end
end
[~,I] = max(obj_fun(:));
[row,col] = ind2sub(size(obj_fun),I); % 2D minimizer
xopt_numerical = xgrid(1,col);
yopt_numerical = ygrid(row,1);
ropt_numerical = obj_fun(row,col);

% *********************************************
%% Compute the optimal solution to the EE maximization equivalent FPP (1) by using Algorithm 1 of the article:
%
%Andrea Pizzo, Alessio Zappone and Luca Sanguinetti, "Solving Energy Efficiency Problems
%through Polynomial Optimization Theory," IEEE Signal Processing Letters, Submitted to.
%
% *********************************************

% *********************************************
% Build the equivalent model
% *********************************************
xstar_poly_vec = NaN*ones(1,maxiter);
ystar_poly_vec = NaN*ones(1,maxiter);
r_poly_vec = NaN*ones(1,maxiter);
lambda_poly_vec = NaN*ones(1,maxiter+1);
lambda_poly = 0; lambda_poly_vec(1) = lambda_poly;
% define optimization varaibles and settings
ops = sdpsettings('solver','sdpt3','sdpt3.gaptol',1e-10,'sdpt3.inftol',1e-10,'sdpt3.steptol',1e-12);

% define the optimization variables
sdpvar x y
n = 2;
% numerator objective
f = (f10bar*x + f11bar*x*y + f20bar*x^2 + f21bar*x^2*y + f30bar*x^3);
% denominator objective
g = (g00 + g10bar*x + g01bar*y + g20bar*x^2 + g02bar*y^2 + g11bar*x*y + g21bar*x^2*y + g12bar*x*y^2 + g30bar*x^3);
%

% *********************************************
% Run Dinkelbach
% ********************************************
for iter=1:maxiter
    
    %Compute the objective function of Problem (3) related to (22)
    p = f-lambda_poly*g;
    %Compute the constraint of (22)
    F = [h00+h10bar*x+h01bar*y>=0, q00+q10bar*x+q01bar*y+q20bar*x^2+q11bar*x*y>=0, x>=1/sK, y>=1/sM, x<=Kmax/sK, y<=Mmax/sM];
    % solve the dual optimization
    solvemoment(F,-p,ops,SOSorderRelaxation);
    value(p)
    solvemoment([F, p<=value(p)],-p,ops,SOSorderRelaxation);
    value(p)
    [sol,~,info] = solvemoment([F, p<=value(p)],-p,ops,SOSorderRelaxation);
    value(p)
    % Lower bound
    pstar_mom = relaxdouble(p);
    % Rank-1 moment matrix from which solution is extracted
    % Simple case, rank-1
    M_mom = info.moment{end};
    S_mom = diag(eig(info.moment{end}));
    % the n first monomials are linear
    sdisplay(info.monomials)';
    assign(info.monomials(2:2+n-1),info.moment{end}(1,(2:2+n-1))');
    
    % return optimizers
    ystar_mom_char = sdisplay(info.monomials(2:2+n-1)');
    ystar_mom = double(info.monomials(2:2+n-1)');
    xstar = ystar_mom(1,1);
    ystar = ystar_mom(1,2);
    check(F)
    
    % Dinkelbach's update rule
    lambda_poly = (f10bar*xstar + f11bar*xstar*ystar + f20bar*xstar^2 + f21bar*xstar^2*ystar + f30bar*xstar^3)...
        /(g00 + g10bar*xstar + g01bar*ystar + g20bar*xstar^2 + g02bar*ystar^2 + g11bar*xstar*ystar + g21bar*xstar^2*ystar + g12bar*xstar*ystar^2 + g30bar*xstar^3);
    
    % Plot obtained minimizer
    rstar_poly = lambda_poly/1e6;
    
    % save optimum
    xstar_poly_vec(iter) = xstar*sK;
    ystar_poly_vec(iter) = ystar*sM;
    lambda_poly_vec(iter+1) = lambda_poly;
    r_poly_vec(iter) = rstar_poly;
        
end


end
