function [PopQs, popstructure] = CR_DEBIPM(timesteps,MatrixSize)

% Please cite:

% To create figure 3:
% [PopQs, popstructure] = CR_DEBIPM(1500,200);

% input: 
% timesteps: number of time steps
% MatrixSize: number of size bins of matrix approximation of DEBIPM

% IPM parameters
Lmin = 150/1000; % minimum body length in mm (used to define limits of matrix approximation)
Lb = 167/1000; % length at birth in mm 
Lm = 1008/1000; % maximum body length in mm
E_Ystdev=0.3; % sigma(Y)

% DEB parameters 
kappa = 0.082; % fraction energy allocation to respiration (as opposed to reproduction)
mu = 0.03; % mortality rate
beta = 151; % to calculate von Bertalanffy growth rate r_b
alpha = -137.8; % to calculate von Bertalanffy growth rate r_b
Rm = 32; % maximum reproduction rate at L_m
d_Lp = 0.539; % proportionality constant for Lp

% Resource parameters
Xmax = 4; % max resource abundance in absence of consumers
Imax = 0.1; % max ingestion rate; default value 0.1 
rho = 0.1; % resource replenishment rate; default value 0.1
K = 5; % half-saturation constant; default value 5

% starting values
X = Xmax;
time = 1;
tick = 1;
ntv = 0.01*ones(MatrixSize, 1); % start population vector
ntv_dead = zeros(1,MatrixSize); % population vector when consumers go extinct

        while time < timesteps, time/timesteps 
            
            E_Y = X/(K+X) ; % functional response

            % calculate IPM and pop quantities
                if sum(ntv)>0, % only calculate kernel if C is alive
                % update DEB pars          
                Linf = E_Y*Lm;
                Lp = max(0.314,d_Lp*Linf); % minimum size for reproduction of 0.314 mm
                rb = 1/(beta + alpha*Linf); % von Bertellanfy growth rate
                % kernel                        
                [S, R, G, D, meshpts] = BigMatrix(MatrixSize,Lmin,Lm,Lp,Rm,Lm,E_Y,E_Ystdev,rb,kappa,Lb,mu);
                kernelDEB = G*S + D*R; 
                % lambda and ssd
                [W,d] = eig(kernelDEB); lambda1 = diag(d); imax = find(lambda1==max(lambda1));
                lambda = max(lambda1(imax)); % population growth rate
                w1=W(:,imax); ssd = w1/sum(w1); ssd = ssd'; % stable stage distribution
                % LRS 
                [m,p]=size(kernelDEB); TS = G*S; Fmat = D*R; 
                R0mat = Fmat*inv((eye(m)-TS)); meanLRS = max(eig(R0mat)); 
                % generation time
                GT = log(meanLRS)/log(lambda);
                %  density of juvs and adults:
                index_adults =  meshpts>=Lp; adultdens = sum(ntv.*index_adults');
                index_juvs =  meshpts<Lp; juvdens = sum(ntv.*index_juvs');
                else disp('C dead'); ntv = ntv_dead'; lambda = 0; meanLRS = 0; GT = 0; Linf = NaN; Lp = NaN; adultdens = NaN; juvdens = NaN; 
                end
            
       popstructure(:,tick) = ntv; 
       PopQstemp = [time,E_Y,sum(ntv),Linf,Lp, adultdens,juvdens, lambda,meanLRS,GT,X];
       PopQs(tick,:) = PopQstemp;
            
       % resource growth
       X_growth = (1-rho)*X + rho*Xmax;
       % resource consumption
       X_consumption = Imax*E_Y*sum(ntv'.*meshpts.^2);
            
       % update population, resource density, time step 
       X = X_growth - X_consumption; 
       X = max(X,0); % to remove negative values
       ntv = kernelDEB*ntv;
       time = time + 1;
       tick = tick + 1;
       end
        
%% Figure
figure
subplot(4,2,1); plot(PopQs(:,1),PopQs(:,11),'k-','LineWidth',2); box off; axis([0 timesteps 0 1.2*max(PopQs(:,11))]); 
ylabel('Resource density','FontSize',14); xlabel('Time'); 
subplot(4,2,3); plot(PopQs(:,1),PopQs(:,2),'k-','LineWidth',2); box off; axis([0 timesteps 0 1.2*max(PopQs(:,2))]); 
ylabel('Feeding level','FontSize',14); xlabel('Time');  
subplot(4,2,5); plot(PopQs(:,1),PopQs(:,7),'k','LineWidth',2); axis([0 timesteps 0 1.2*max(PopQs(:,7))]); box off
ylabel('Juv density','FontSize',14);  xlabel('Time');  set(gca,'FontSize',12); 
subplot(4,2,7); plot(PopQs(:,1),PopQs(:,6),'k','LineWidth',2); axis([0 timesteps 0 1.2*max(PopQs(:,6))]); box off
ylabel('Adult density','FontSize',14);  xlabel('Time');  set(gca,'FontSize',12); 
        