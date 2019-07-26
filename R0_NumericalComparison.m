%---------------------------------------------------------------------------
%% R0_NumericalComparison.m
%% Author: Leah Darwin, Arizona State University 
%% Contact: ljdarwin@asu.edu
%% Description: MATLAB code to produce "Numerical Comparison of R_0 Values"
%% used in paper "The Effects of Pharmacological Chaperone and Interferon Treatment on Prion Proliferation in the Brain"
%% shows how different analytical values for R_0 vary numerically.  Produces a single figure. 
%% Rights Statement:  All rights reserved MTBI 2019. 
%% Last Edits Made: July 24, 2019
%---------------------------------------------------------------------------

%%house keeping
clear all
close all
hold off

%%typeset LaTeX interperter
set(0,'defaulttextinterpreter','latex')

%%pre-declaration of parameter values 
%---------------------------------------------------------------------------
Lambda = 2400; 
I = 8.94849*10^6/5;
mu_P = 0.047;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_S = 4; mu_A = 62.0352; mu_R=4; mu_I = k*mu_P;
beta_S = 2.29*10^-3; beta_R = beta_S; 
b = 0.0314;
n = 3;
sigma = 1;
alpha = 5.1408*10^-2; 
%---------------------------------------------------------------------------

%%declare vectors of treatment samples
D = 0:(500000/7)/1000:500000/7;

%%S PFE calculation
S_PFE = 1/2.*(Lambda./mu_S-D./(mu_A+mu_S)-mu_A./alpha)+sqrt((Lambda*mu_A)./(alpha.*mu_S)+1/4*(D./(mu_A+mu_S)+mu_A./alpha-Lambda./mu_S).^2);

%%Heuristic R0 calculation
R0_heuristic=(sqrt(beta_S.*b.*S_PFE+b.^2/4)-b./2)./((mu_P+mu_I)+(b.*(n-1)));

%%Next Generation Matrix R0 calculation
R0_nextGen=(b.*(beta_S.*S_PFE-n.*(n-1).*b))./((mu_P+mu_I).*((mu_P+mu_I)+(2.*n-1).*b));

%%Square Root of Next Generation Matrix R0 calculation
nextGen_sqrt = sqrt(R0_nextGen);
subplot(2,1,1)

%%plot all 3 R0 values with variation in pharamcological chaperone dosage
plot(D,R0_heuristic,D,R0_nextGen,':', D,nextGen_sqrt,'--','LineWidth',3)
title('\textbf{Numerical Comparisons for $R_0$}','Interpreter','latex','FontSize',30);
legend({'$R_0^H$','$R_0^{NG}$','$\sqrt{R_0^{NG}}$'},'Interpreter','latex','FontSize',20)
xlabel('Dose of Interferons (nM/day)','Interpreter','latex','FontSize',20)
ylabel('Value of $R_0$')
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
axis tight
ax = gca;
ax.FontSize = 20; 

%%update parameters
I = 0:(8.94849*10^6/2)/1000:8.94849*10^6/2;
k=-log(1./(1+(I./(1000.*13989.5)).^1.44313))./(2.*mu_P);
mu_I = k.*mu_P;
D = 500000/7;

%%S PFE calculation
S_PFE = 1/2.*(Lambda./mu_S-D./(mu_A+mu_S)-mu_A./alpha)+sqrt((Lambda*mu_A)./(alpha.*mu_S)+1/4*(D./(mu_A+mu_S)+mu_A./alpha-Lambda./mu_S).^2);

%%Heuristic R0 calculation
R0_heuristic=(sqrt(beta_S.*b.*S_PFE+b.^2/4)-b./2)./((mu_P+mu_I)+(b.*(n-1)));

%%Next Generation Matrix R0 calculation
R0_nextGen=(b.*(beta_S.*S_PFE-n.*(n-1).*b))./((mu_P+mu_I).*((mu_P+mu_I)+(2.*n-1).*b));

%%Square Root of Next Generation Matrix R0 calculation
nextGen_sqrt = sqrt(R0_nextGen);

%%plot all 3 R0 values with variation in interferon dosage
subplot(2,1,2)
plot(mu_I, R0_heuristic, mu_I, R0_nextGen,':', mu_I,nextGen_sqrt,'--','LineWidth',3)
legend({'$R_0^H$','$R_0^{NG}$','$\sqrt{R_0^{NG}}$'},'Interpreter','latex','FontSize',20)
xlabel('Dose of PCs (nM/day)','Interpreter','latex','FontSize',20)
ylabel('Value of $R_0$')
axis tight

%%format final graph
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
axis tight
ax = gca;
ax.FontSize = 20; 


