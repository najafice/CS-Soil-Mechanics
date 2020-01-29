%% Matlab Code for Cam-Clay Model (Drained and Undrained)
%% Author: Abolfazl Najafi, Iran University of S & T.
% Under supervision of Dr. Alireza Saeidi-Azizkandi, Ph.D

clc;
clear all;
close all;

%% Display Product Name, Author and License Information
    disp(' ');
    disp('Matlab code for Simulation of CamClay model');
    disp('Abolfazl Najafi, Iran University of S & T.');
    
%% Input Parameters
    disp(' ');
    disp('Input Parameters for Cam-Clay:');
    
    cp=input('Enter the inital Consolidation pressure (kPa) (eg., 100 kPa)  = ');
    dp=input('Enter the value of Deviatoric Stress (kPa)    (eg., 300 kPa)  = ');
    M=input( 'Enter the value of Critical Friction Angle M  (eg., 0.95)     = ');
    N=input( 'Enter the value of N                          (eg., 2.5)      = ');
    l=input( 'Enter the value of Lamda                      (eg., 0.093)    = ');
    k=input( 'Enter the value of Kappa                      (eg., 0.035)    = ');
    
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

%% Computation of Other Parameters (V0 and Gamma)
    V0=N-(l*log(cp));               % Initial Specific Volume
    G=N-(l-k);                      % Value of Gamma constant

%% Strain Increament and Strain Matrix Definition
    iter=2;

%% Block Memory allocation
    u=zeros(iter,1);    % Pore Water Pressure
    p=zeros(iter,1);    % Mean Effective Stress
    q=zeros(iter,1);    % Deviatoric Stress
    epsV=zeros(iter,3); % Volumetic Strain
    epsS=zeros(iter,1); % Shear Strain
    
%% Initialize   
     S=[cp,cp,cp;dp,cp,cp];       % Stress 
     
     p(1)=(S(1,1)+2*S(1,3))/3;                
     q(1)=(S(1,1)-S(1,3));
     
     p(2)=(S(2,1)+2*S(2,3))/3;                
     q(2)=(S(2,1)-S(2,3));
     
%% CamClay Iteration Uni-Loop Iteration for OC/NC & Inside/Outside Yield
     for i=2:iter

       %Stress and Strain Updates
       if analysis==1 %Triaxial Drained
           deltap=p(i)-p(i-1);
           deltaq=q(i)-q(i-1);
           
           epsV(i,1)=(((l-k)/(M*V0*cp))*((M-(q(1)/p(1)))*deltap+deltaq)); % Plastic Volumetric Strain
           epsV(i,2)=(k*deltap/(V0*p(1))); % Elastic Volumetric Strain
           epsV(i,3)=epsV(i,1)+epsV(i,2); % Total Value of Volumetric Strain
                       
           epsS(i,1)=((1/(M-(q(1)/p(1))))*epsV(i,1)); % Plastic Shear Strain
           
           disp(['Volumetric Strain is: ', num2str(epsV(i,3)*100)]);
           disp(['Shear Strain is: ', num2str(epsS(i,1)*100)]);
       elseif analysis==2 %Triaxial Undrained
           syms pp
           eqn = q(2) == (M*pp/(l-k))*(G+l-k-V0-l.*log(pp));
           pprime=real(vpasolve(eqn,pp, [10 p(2)]));
           
           deltapp=pprime-p(i-1);
           deltap=p(i)-p(i-1);
           deltaq=q(i)-q(i-1);
           
           u(i)=deltap-deltapp;
           
           epsV(i,1)=-(k*deltapp)/(V0*cp); % Plastic Volumetric Strain
           epsS(i,1)=((1/(M-(q(1)/p(1))))*epsV(i,1)); % Plastic Shear Strain
           
           disp(['Volumetric Strain is: ', num2str(epsV(i,1)*100)]);
           disp(['Shear Strain is: ', num2str(epsS(i,1)*100)]);
           disp(['Pore Pressure is: ', num2str(u(i))]);
       end    
     end