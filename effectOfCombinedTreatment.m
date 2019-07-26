%---------------------------------------------------------------------------
%% effectOfCombinedTreatment.m
%% Author: Leah Darwin, Arizona State University 
%% Contact: ljdarwin@asu.edu
%% Description: MATLAB code to produce "Effects of Treatment on Prion Proliferation"
%% used in paper "The Effects of Pharmacological Chaperone and Interferon Treatment on Prion Proliferation in the Brain"
%% shows how prion dynamics are affected by the introduction of combined treatments.  Produces a single figure. 
%% Rights Statement:  All rights reserved MTBI 2019. 
%% Last Edits Made: July 24, 2019
%---------------------------------------------------------------------------

%%house keeping
clear all
close all
hold off

%%%%forward declaration of parameter values (subject to change)
Lambda = 2400; 
I = 0;
mu_P = 0.047;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_S = 4; mu_A =  62.0352; mu_R=4; mu_I = k*mu_P*0;
beta_S = 2.29*10^-3 ; beta_R = beta_S ; 
b = 0.0314;
n = 3;
sigma = 1;
alpha =  5.1408*10^-2; 
D = 0;

%%declare time vector (~1 year)
time = 0:400;

%%NO TREATMENT SECTION
%%**************************************************************************
%%declare/update system of differential equations
%%let the vector x = [ S, T, P, Z, R, A ]
f = @(t,x) [ Lambda-mu_S*x(1)-alpha*x(6)*x(1)+mu_A*x(2)-beta_S*x(1)*x(3);
            -mu_S*x(2)+alpha*x(6)*x(1)-mu_A*x(2);
            -(mu_P+mu_I)*x(3)+b*x(4)-(2*n-1)*b*x(3);
            beta_S*x(1)*x(3)+beta_R*x(5)*x(3)-(mu_P+mu_I)*x(4)-n*(n-1)*b*x(3);
            -mu_R*x(5)-mu_I*sigma*x(5)-beta_R*x(5)*x(3)+n*(n-1)*b*x(3);
            -mu_A*x(6)-alpha*x(6)*x(1)+D];

%%solve the system of differential equations        
[t,xa] = ode45(f,time,[750,0,3000,3*30,0,0]);

%%plot P state variable against time
plot(time,(xa(:,3))','LineWidth',3);
hold on 
%%**************************************************************************

%%JUST PCs SECTION
%%**************************************************************************
%%update parameter values
I = 0;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_I = k*mu_P*0;
D = (500000/7);

%%declare/update system of differential equations
%%let the vector x = [ S, T, P, Z, R, A ]
f = @(t,x) [ Lambda-mu_S*x(1)-alpha*x(6)*x(1)+mu_A*x(2)-beta_S*x(1)*x(3);
            -mu_S*x(2)+alpha*x(6)*x(1)-mu_A*x(2);
            -(mu_P+mu_I)*x(3)+b*x(4)-(2*n-1)*b*x(3);
            beta_S*x(1)*x(3)+beta_R*x(5)*x(3)-(mu_P+mu_I)*x(4)-n*(n-1)*b*x(3);
            -mu_R*x(5)-mu_I*sigma*x(5)-beta_R*x(5)*x(3)+n*(n-1)*b*x(3);
            -mu_A*x(6)-alpha*x(6)*x(1)+D];

%%solve the system of differential equations  
[t,xa] = ode45(f,time,[750,0,3000,3*30,0,0]);

%%plot P state variable against time
plot(time,(xa(:,3)),'--','LineWidth',3);
hold on 
%%**************************************************************************

%%JUST INTERFERONS
%%**************************************************************************
%%update parameter values
I = ((5*10^6)/2)/2;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_I = k*mu_P;
D = 0;

%%declare/update system of differential equations
%%let the vector x = [ S, T, P, Z, R, A ]
f = @(t,x) [ Lambda-mu_S*x(1)-alpha*x(6)*x(1)+mu_A*x(2)-beta_S*x(1)*x(3);
            -mu_S*x(2)+alpha*x(6)*x(1)-mu_A*x(2);
            -(mu_P+mu_I)*x(3)+b*x(4)-(2*n-1)*b*x(3);
            beta_S*x(1)*x(3)+beta_R*x(5)*x(3)-(mu_P+mu_I)*x(4)-n*(n-1)*b*x(3);
            -mu_R*x(5)-mu_I*sigma*x(5)-beta_R*x(5)*x(3)+n*(n-1)*b*x(3);
            -mu_A*x(6)-alpha*x(6)*x(1)+D];

%%solve the system of differential equations  
[t,xa] = ode45(f,time,[750,0,3000,3*30,0,0]);

%%plot P state variable against time
plot(time,(xa(:,3))',':','LineWidth',3);
hold on 
%%**************************************************************************

%%COMBINED TREATMENT
%%**************************************************************************
%%update parameter values
I = ((5*10^6)/2)/2;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_I = k*mu_P;
D = (500000/7);

%%declare/update system of differential equations
%%let the vector x = [ S, T, P, Z, R, A ]
f = @(t,x) [ Lambda-mu_S*x(1)-alpha*x(6)*x(1)+mu_A*x(2)-beta_S*x(1)*x(3);
            -mu_S*x(2)+alpha*x(6)*x(1)-mu_A*x(2);
            -(mu_P+mu_I)*x(3)+b*x(4)-(2*n-1)*b*x(3);
            beta_S*x(1)*x(3)+beta_R*x(5)*x(3)-(mu_P+mu_I)*x(4)-n*(n-1)*b*x(3);
            -mu_R*x(5)-mu_I*sigma*x(5)-beta_R*x(5)*x(3)+n*(n-1)*b*x(3);
            -mu_A*x(6)-alpha*x(6)*x(1)+D];

%%solve the system of differential equations  
[t,xa] = ode45(f,time,[750,0,3000,3*30,0,0]);

%%plot P state variable against time
plot(time,(xa(:,3))','-.','LineWidth',3);
hold on 
%%**************************************************************************

%%format final plot
%---------------------------------------------------------------------------
xlabel('Time (days)','Interpreter','latex','FontSize',30)
ylabel('Prion Concentration (nM)','Interpreter','latex','FontSize',30)
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis tight
ax = gca;
ax.FontSize = 20; 
title('\textbf{Effects of Treatment on Prion Proliferation}','Interpreter','latex','FontSize',40)
lgd = legend({'No Treatment','Pharmacological Chaperones','Interferons','Combined'},'Interpreter','latex','Location','northwest');
lgd.FontSize = 24;
%---------------------------------------------------------------------------




