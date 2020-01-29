%% Matlab Code for Roscoe-Worslef Surfaces
%% Author: Abolfazl Najafi, Iran University of S & T.
% Under supervision of Dr. Alireza Saeidi-Azizkandi, Ph.D

clc;
clear all;
close all;

%% Display Product Name, Author and License Information
    disp(' ');
    disp('Matlab code for plotting Roscoe & Worslef surfaces.');
    disp('Abolfazl Najafi, Iran University of S & T.');
    
%% Input Parameters
    disp(' ');
    disp('Input Parameters for Roscoe-Worslef:');
    
     M=input( 'Enter the value of Critical Friction Angle M  (eg., 1.02)     = ');
     N=input( 'Enter the value of N                          (eg., 3.32)      = ');
     g=input( 'Enter the value of Gamma                      (eg., 3.17)      = ');
     l=input( 'Enter the value of Lamda                      (eg., 0.2)      = ');
     k=input( 'Enter the value of Kappa                      (eg., 0.05)     = ');
     h=input( 'Enter the value of h                          (eg., 0.675)    = ');
    
%    M=1.02; N=3.32; g=3.17; Landa=0.2; k=0.05; h=0.675;
    
    n=50;
%% Analysis and Plotting section

    v1=linspace(1,N,n);

    pnc=exp((N-v1)./l);
    pcs=exp((g-v1)./l);
    ptt=((M-h)/(3-h))*exp((g-v1)/l);

    for i=1:length(pcs)
        pr(i,:)=linspace(pnc(i),pcs(i),n);
    end

    for i=1:length(pcs)
        pv(i,:)=linspace(ptt(i),pcs(i),n);
    end

    for i=1:length(pcs)
        pt(i,:)=linspace(ptt(i),0,n);
    end

    v=repmat(v1,n,1)';
    
    qr=((M.*pr)./(l-k)).*(g+l-k-v-l.*log(pr));
    qv=(M-h)*exp((g-v)/l)+h*pv;
    qt=3.*pt;
       
    b=surf(pr,v,qr);
    hold on
    c=surf(pv,v,qv);
    d=surf(pt,v,qt);
    plot3(pcs',v(:,1),qv(:,n))
    plot3(ptt',v(:,1),qt(:,1))
    plot3(pnc',v(:,1),qr(:,1).*0)
    % hold off
    grid on
    hAxis = gca;
    hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
    hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
    hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
    hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis

    hAxis.ClippingStyle = 'rectangle';
    b.MeshStyle = 'row';
    b.FaceColor = 'none';

    c.MeshStyle = 'row';
    c.FaceColor =  'none';

    d.MeshStyle = 'row';
    d.FaceColor =  'none';

    xlim([0 inf])
    ylim([0.9 N])
    zlim([0 inf])
    view([-130 20])