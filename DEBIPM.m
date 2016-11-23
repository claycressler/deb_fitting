function [E_Y, PopQs] = DEBIPM(E_Ymin,E_Ymax,step,MatrixSize)

% Please cite:

% To create figure 2:
% [E_Y, PopQs] = DEBIPM(0.65,0.95,15,200);

% input: 
% E_Ymin: minimum feeding level
% E_Ymax: maximum feeding level
% step: step size between E_Ymin and E_Ymax
% MatrixSize: number of size bins of matrix approximation of DEBIPM

Lmin = 150/1000; % minimum body length in mm (used to define limits of matrix approximation)
Lb = 167/1000; % length at birth in mm 
Lm = 1008/1000; % maximum body length in mm 
E_Ystdev=0.1; % sigma(Y)

% DEB parameters 
kappa = 0.082; % fraction energy allocation to respiration (as opposed to reproduction)
mu = 0.03; % mortality rate
beta = 151; % to calculate von Bertalanffy growth rate r_b
alpha = -137.8; % to calculate von Bertalanffy growth rate r_b
Rm = 32; % maximum reproduction rate at L_m

E_Y=linspace(E_Ymin,E_Ymax,step); % Feeding level
       
        for n=1:length(E_Y); 
            
            Linf = E_Y(n)/E_Ymax*Lm; % ultimate length  
            Lp = 0.539*Linf; % length at puberty
            rb = 1/(beta + alpha*Linf); % von Bertellanfy growth rate
                        
            [S, R, G, D, meshpts] = BigMatrix(MatrixSize,Lmin,Linf,Lp,Rm,Lm,E_Y(n),E_Ystdev,rb,kappa,Lb,mu);
            kernel = G*S + D*R; 
            % calculate lambda, ssd, and rv
            [W,d] = eig(kernel); lambda1 = diag(d); imax = find(lambda1==max(lambda1)); 
            V=conj(inv(W)); lambda = lambda1(imax); % population growth rate
            w1=W(:,imax); v1 = real(V(imax,:))';
            ssd = w1/sum(w1); ssd = ssd'; % stable stage distribution

            % calculation of LRS
            [m,p]=size(kernel); TS = G*S; Fmat = D*R; 
            R0mat = Fmat*inv((eye(m)-TS)); meanLRS = max(eig(R0mat)); 

            %  generation time
            GT = log(meanLRS)/log(lambda);

            PopQstemp = [E_Y(n),Lm,Lp,Rm,lambda,meanLRS,GT];
            PopQs(n,:) = PopQstemp;
        end
        
        subplot(3,1,1); plot(E_Y,PopQs(:,5),'k-','LineWidth',2);  axis square;  box off
        ylabel('Population growth rate','FontSize',14); xlabel(' ');  
        subplot(3,1,2); plot(E_Y,PopQs(:,6),'k-','LineWidth',2);   axis square;  box off
        ylabel('R_0','FontSize',14);  xlabel(' ');  set(gca,'FontSize',12); 
        subplot(3,1,3); plot(E_Y,PopQs(:,7),'k-','LineWidth',2); axis square;  box off; hold on
        ylabel('Generation time','FontSize',14);  xlabel('Feeding level');  set(gca,'FontSize',12); 