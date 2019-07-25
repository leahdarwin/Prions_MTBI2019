%---------------------------------------------------------------------------
%% effectOfPCDosage.m
%% Author: Leah Darwin, Arizona State University 
%% Contact: ljdarwin@asu.edu
%% Description: MATLAB code to produce "Effect of Pharmacological Chaperon Dosage on Prion Dynamics"
%% used in paper "The Effects of Pharmacological Chaperone and Interferon Treatment on Prion Proliferation in the Brain"
%% shows how prion dynamics are affected by the introduction of Interferon treatment.  Produces a single figure. 
%% Rights Statement:  All rights reserved MTBI 2019. 
%% Last Edits Made: July 24, 2019
%---------------------------------------------------------------------------

%%house keeping
clear all
close all
hold off

%%typeset LaTeX interperter
set(0,'defaulttextinterpreter','latex')

%%forward declaration of parameter values (subject to change)
%---------------------------------------------------------------------------
Lambda = 2400; 
I = 0;
mu_P = 0.047;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_S = 4; mu_A = 62.0352; mu_R=4; mu_I = k*mu_P*0;
beta_S = 2.29*10^-3 ; beta_R = beta_S ; 
b = 0.0314;
n = 3;
sigma = 1;
alpha = 5.1408*10^-2; 
D = 0;
%---------------------------------------------------------------------------

%%declaration of time vector (2 years)
time = 0:365*2;

%%declare vector of interferon dosages to sample from
D_samples = [0:(500000/7)/4:500000/7 500000/7+(500000/7)/4:(500000/7)/4:500000/2];

%%preallocation of P state variable equilbrium values
final_P = zeros(1,length(D_samples));

%%code snipit to create color gradient for plot lines
%---------------------------------------------------------------------------
%%create vector to store color for every plot line that will be created
lengthColor = length(D_samples);
%%intialize starting color of gradient
color1 = [0, 0, 1];
%%intialize ending color of gradients
color2 = [126, 206, 253]/255;
%%interpolate all colors inbetween start and end color
colors_p = [linspace(color1(1),color2(1),lengthColor)', linspace(color1(2),color2(2),lengthColor)', linspace(color1(3),color2(3),lengthColor)'];
%---------------------------------------------------------------------------

%%Solve the system of differentials for every value of D sampled
for i=1:length(D_samples)
    
%%update paramter values for each D        
D = D_samples(i);

%%declare/update system of differential equations
%%let the vector x = [ S, T, P, Z, R, A ]
f = @(t,x) [ Lambda-mu_S*x(1)-alpha*x(6)*x(1)+mu_A*x(2)-beta_S*x(1)*x(3);
            -mu_S*x(2)+alpha*x(6)*x(1)-mu_A*x(2);
            -(mu_P+mu_I)*x(3)+b*x(4)-(2*n-1)*b*x(3);
            beta_S*x(1)*x(3)+beta_R*x(5)*x(3)-(mu_P+mu_I)*x(4)-n*(n-1)*b*x(3);
            -mu_R*x(5)-mu_I*sigma*x(5)-beta_R*x(5)*x(3)+n*(n-1)*b*x(3);
            -mu_A*x(6)-alpha*x(6)*x(1)+D];

%%solve the system of differential equations
[t,xa] = ode45(f,time,[750,0,3,3*30,0,0]);

%%update P state variable equilbrium solution vector
final_P(i)= xa(length(xa),3);

%%plot P against time for given value of D
subplot(2,1,1);
%%if D is below or at toxic threshold, plot with a solid line
if D <= 500000/7
    plot(time,((xa(:,3)))','color',colors_p(i,:),'LineWidth',3);
    hold on
%%if D is above the toxic threshold, plot with a dotted line
else
    plot(time,((xa(:,3)))','--','color',colors_p(i,:),'LineWidth',3);
    hold on
end
end

%%format D vs Time graph
%---------------------------------------------------------------------------
xlabel('Time (days)','Interpreter','latex','FontSize',20)
ylabel('Prion Concentration (nM)','Interpreter','latex','FontSize',20)
title('\textbf{Effect of PCs on Prion Dynamics}','Interpreter','latex','FontSize',60);
ax = gca;
ax.FontSize = 16; 
%%declaration/formatting of color bar
colormap(colors_p);
c = colorbar('Ticks',[0 1],'TickLabels',[D_samples(1) D_samples(length(D_samples))],'TickLabelInterpreter','latex','FontSize',16);
c.Label.String = 'PC Dose (nM/day)'
c.Label.Interpreter = 'latex';
axis tight
%---------------------------------------------------------------------------

%plot P state variable equilbrium solution vector against time
%---------------------------------------------------------------------------

%%create subvectors for toxic and nontoxic doses of D 
D_samples_nonToxic = [];

%%divide D vector non-toxic vector
for i=1:length(D_samples)
    if D_samples(i) <= 500000/7
        D_samples_nonToxic = [D_samples_nonToxic D_samples(i)];
    end
end

%%plot all doses with a dotted line
subplot(2,1,2);
plot(D_samples,final_P,'--','color',[0,0,1],'LineWidth',3);
hold on

%plot an overlaid solid line for nontoxic doses 
subplot(2,1,2);
plot(D_samples_nonToxic,final_P(1,1:length(D_samples_nonToxic)),'color',[0,0,1],'LineWidth',3);
hold on

%%format P_final vs Time graph
%---------------------------------------------------------------------------
xlabel('PC Dose (nM/day)','Interpreter','latex','FontSize',20)
ylabel('Final Prion Concentration (nM)','Interpreter','latex','FontSize',20)
axis tight
ax = gca;
ax.FontSize = 20;
%---------------------------------------------------------------------------

%%plot toxic dosage line
%---------------------------------------------------------------------------
%%declare y range vector
y = 0:100:5000;
%%make toxic x values to match y length
toxic_vec = ones(length(y))*500000/7;
plot(toxic_vec,y,'color',[1 0 0],'LineWidth',3);
%---------------------------------------------------------------------------