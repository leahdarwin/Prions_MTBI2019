%---------------------------------------------------------------------------
%% effectOfInterferonDosage.m
%% Author: Leah Darwin, Arizona State University 
%% Contact: ljdarwin@asu.edu
%% Description: MATLAB code to produce "Effect of Interferon Dosage on Prion Dynamics"
%% used in paper "The Effects of Pharmacological Chaperone and Interferon Treatment on Prion Proliferation in the Brain"
%% shows how prion dynamics are affected by the introduction of Interferon treatment.  Produces a single figure. 
%% Rights Statement:  All rights reserved MTBI 2019. 
%---------------------------------------------------------------------------

%%house keeping
clear all
close all
hold off

%typeset LaTeX interperter
set(0,'defaulttextinterpreter','latex')


%%forward declaration of parameter values (subject to change)
%---------------------------------------------------------------------------
Lambda = 2400; 
I = 0;
mu_P = 0.047;
k=-log(1/(1+((I*2)/(1000*13989.5))^1.44313))/(2*mu_P);
mu_S = 4; mu_A = 62.0352; mu_R=4; mu_I = k*mu_P;
beta_S = 2.29*10^-3 ; beta_R = beta_S ; 
b = 0.0314;
n = 3;
sigma = 1;
alpha = 5.1408*10^-2; 
D = 0;
%---------------------------------------------------------------------------

%%declaration of time vector (1 year)
time = 0:365;

%%declare vector of interferon dosages to sample from
I_samples = 0:(8.94849*10^6/2)/10:8.94849*10^6/2;

%%preallocation of P state variable equilbrium values
final_P = zeros(1,length(I_samples));

%%code snipit to create color gradient for plot lines
%---------------------------------------------------------------------------
%%create vector to store color for every plot line that will be created
lengthColor = length(I_samples);
%%intialize starting color of gradient
color1 = [1, 0, 0];
%%intialize ending color of gradients
color2 = [255, 255, 0]/255;
%%interpolate all colors inbetween start and end color
colors_p = [linspace(color1(1),color2(1),lengthColor)', linspace(color1(2),color2(2),lengthColor)', linspace(color1(3),color2(3),lengthColor)'];
%---------------------------------------------------------------------------

%%Solve the system of differentials for every value of I sampled
for i=1:length(I_samples)
    
%%update paramter values for each I    
I = I_samples(i);
k=-log(1/(1+(2*I/(1000*13989.5))^1.44313))/(2*mu_P);
mu_I = k*mu_P;

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

%%plot P against time for given value of I
subplot(2,1,1);
plot(time,((xa(:,3)))','color',colors_p(i,:),'LineWidth',3);
hold on
end

%%format P vs Time graph
%---------------------------------------------------------------------------
axis tight
xlabel('Time (days)','Interpreter','latex','FontSize',20)
ylabel('Prion Concentration (nM)','Interpreter','latex','FontSize',20)
title('\textbf{Effect of Interferon Dosage on Prion Dynamics}','Interpreter','latex','FontSize',60);
ax = gca;
ax.FontSize = 16; 
%%declaration/formatting of color bar
colormap(colors_p);
c = colorbar('Ticks',[0 1],'TickLabels',[I_samples(1) I_samples(length(I_samples))],'TickLabelInterpreter','latex','FontSize',16);
c.Label.String = 'Interferon Dose (nM/day)'
c.Label.Interpreter = 'latex';
%---------------------------------------------------------------------------

%plot P state variable equilbrium solution vector against time
subplot(2,1,2);
plot(I_samples,final_P,'color',[1,0,0],'LineWidth',3);

%%format P_final vs Time graph
%---------------------------------------------------------------------------
xlabel('Interferon Dose (nM/day)','Interpreter','latex','FontSize',20)
ylabel('Final Prion Concentration (nM)','Interpreter','latex','FontSize',20)
axis tight
ax = gca;
ax.FontSize = 20; 
%---------------------------------------------------------------------------


