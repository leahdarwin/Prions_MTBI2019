%%GrRavsIn&Pcs.m
%%Author: Yair Antonio Castillo-Castillo, University of Colima
%%Contact: dekelblack@gmail.com
%%Description: Matlab code to produce "Growth rate vs concentration of
%%treatments used in the paper "The Effects of Pharmacological Chpaeron and
%%Interferon Treatment on Prion Proliferation in the Brain" shows how the
%%growth rate is affected by ony each concentration of Interferons and
%%Pharmacological Treatments. Produce a single figure-
%%Rights Statements: All rights reserved MTBI 2019
%%Last Edits Made: July 24, 2019
%% Intervals
c=[0:1000:5000000]; %Interval goes to 0 to 5000000 

%% FUNCTION OF GROWTH RATE THAT WE NEED
f=-0.1255+0.0886002.*sqrt(-2.06125-0.0000693569.*c+0.00916.*sqrt(957375+(3.4596+0.0000573309.*c).*c)); % Function of Pharmacological Chaperons
ff=0.5.*(0.214515 + log(1./(1 + 1.32737*10.^-10.*c.^1.443))); % Function of Interferons

%% PLOTS
plot(c,f,'k--')                                                          %Plot Growth rate vs PCs
hold on                                                                  %To plot in the same graph
plot(c,ff,'r-')                                                          %Plot growth rate vs Interferons
line([0 5*10^6], [0 0],'Color','black','LineStyle','-','LineWidth',1.5); %Plot a line (x-axes)
line([500000/7 500000/7], [0.2 -0.2]);                                   %Plot a line (Threshold for Interferons) 
 
%% LABELS FOR THE GRAPHS
hold on
xlabel('Concentration (nM/day)','Interpreter','latex','FontSize',25) %Make the x-axes label
ylabel('Growth rate (r)','Interpreter','latex','FontSize',25)        %Make the y-axes label
title('\textbf{Growth rate vs concentration of treatments}','Interpreter','latex','FontSize',40);                     %Title
legend({['Pharmacological Chaperones'],['Interferons']}, 'Interpreter','latex','Location','northeast','FontSize',25); %Legend