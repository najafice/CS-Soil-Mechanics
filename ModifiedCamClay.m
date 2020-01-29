%% Matlab Code for Modified Cam-Clay Model (Drained and Undrained)
%% Author: Abolfazl Najafi, Iran University of S & T.
% Under supervision of Dr. Alireza Saeidi-Azizkandi, Ph.D

clc;
clear all;
close all;

%% Display Product Name, Author and License Information
    disp(' ');
    disp('Matlab code for Simulation of Modified CamClay model');
    disp('Abolfazl Najafi, Iran University of S & T.');
    
%% Input Parameters
    disp(' ');
    disp('Input Parameters for Modified Cam-Clay:');
    
    cp=input('Enter the inital Consolidation pressure (kPa) (eg., 100 kPa)  = ');
    dp=input('Enter the value of Deviatoric Stress (kPa)    (eg., 300 kPa)  = ');
    M=input( 'Enter the value of Critical Friction Angle M  (eg., 0.95)     = ');
    N=input( 'Enter the value of N                          (eg., 2.5)      = ');
    l=input( 'Enter the value of Lamda                      (eg., 0.093)    = ');
    k=input( 'Enter the value of Kappa                      (eg., 0.035)    = ');
    nu=input('Enter the value of poissons ratio             (eg., 0.15)     = ');
    
    analysis = input('Enter the type of Analysis: (1) Triaxial Drained (2) Triaxial Undrained = ');
    
    disp(' ');
    if analysis==1      % Triaxial Drained
      disp('Triaxial Drained Simulation is in progress ...');
    elseif analysis==2 % Triaxial Undrained
	   disp('Triaxial Undrained Simulation is in progress ...');
	 else           
	   disp('Matlab code handles only Triaxial Drained and Undrained simulations!');
	   quit;
    end 
    
%% Computation of Other Parameters (V0 and G)
    V0=N-(l*log(cp));               % Initial Specific Volume
    
%% Strain Increament and Strain Matrix Definition
    iter=2;
    
%% Block Memory allocation
    u=zeros(iter,1);    % Pore Water Pressure
    p=zeros(iter,1);    % Mean Effective Stress
    q=zeros(iter,1);    % Deviatoric Stress
    epsV=zeros(iter,3); % Volumetic Strain
    epsS=zeros(iter,3); % Shear Strain
    
%% Initialize   
     S=[cp,cp,cp;dp,cp,cp];       % Stress 
     
     p(1)=(S(1,1)+2*S(1,3))/3;                
     q(1)=(S(1,1)-S(1,3));
     
     p(2)=(S(2,1)+2*S(2,3))/3;                
     q(2)=(S(2,1)-S(2,3));
     
%% CamClay Iteration Uni-Loop Iteration for OC/NC & Inside/Outside Yield
     for i=2:iter
       KM=V0*p(i)/k;                     % Bulk Modulus
       G=(3*KM*(1-2*nu))/(2*(1+nu));     % Shear Modulus
       %Stress and Strain Updates
       if analysis==1 %Triaxial Drained
           deltap=p(i)-p(i-1);
           deltaq=q(i)-q(i-1);
           eta=q(i)/p(i);
           
           epsV(i,1)=((l-k)/(V0*p(i)*(M^2+eta^2)))*((M^2-eta^2)*deltap + (2*eta)*deltaq); % Plastic Volumetric Strain
           epsV(i,2)=(k/(V0*p(i)))*deltap; % Elastic Volumetric Strain
           epsV(i,3)=epsV(i,1)+epsV(i,2); % Total Value of Volumetric Strain
                       
           epsS(i,1)=((l-k)/(V0*p(i)*(M^2+eta^2)))*(((2*eta)*deltap) + (((4*eta^2)/(M^2-eta^2))*deltaq)); % Plastic Shear Strain
           epsS(i,2)=(1/(3*G))*deltaq; % Elastic Shear Strain
           epsS(i,3)=epsS(i,1)+epsS(i,2);
           
           disp(['Volumetric Strain is: ', num2str(epsV(i,3)*100)]);
           disp(['Shear Strain is: ', num2str(epsS(i,3)*100)]);
       elseif analysis==2 %Triaxial Undrained
           eta=q(i)/p(i);
           
           syms pp
           eqn = q(i)^2 + (M^2)*pp^2 == (M^2)*pp*p(i-1)^2;
           pprime=max(double(solve(eqn)));
           %pprime=real(vpasolve(eqn,pp, [10 p(i)]));
           
           deltapp=pprime-p(i-1);
           deltap=p(i)-p(i-1);
           deltaq=q(i)-q(i-1);
           
           u(i)=deltap-deltapp;
           epsV(i,1)=((l-k)/((nu*pprime)*(M^2+eta^2))) * (((M^2-eta^2)*deltap) + ((2*eta)*deltaq)); % Plastic Volumetric Strain
           epsS(i,1)=((l-k)/((nu*pprime)*(M^2+eta^2))) * ( ((2*eta)*deltap) + ((4*eta^2)/(M^2-eta^2))*deltaq ); % Plastic Shear Strain
           
           epsV(i,2)=(k/(V0*p(i)))*deltap; % Elastic Volumetric Strain
           epsS(i,2)=(1/(3*G))*deltaq; % Elastic Shear Strain
           
           epsV(i,3)=epsV(i,1)+epsV(i,2); % Total Value of Volumetric Strain
           epsS(i,3)=epsS(i,1)+epsS(i,2); % Total Value of Shear Strain
           
           disp(['Volumetric Strain is: ', num2str(epsV(i,3)*100)]);
           disp(['Shear Strain is: ', num2str(epsS(i,3)*100)]);
           disp(['Pore Pressure is: ', num2str(u(i))]);
       end    
     end