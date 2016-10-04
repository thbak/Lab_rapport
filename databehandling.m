clc
clear all
close all

% importerer al data
data{1} = importdata('100Hz.txt') ;
data{2} = importdata('200Hz.txt') ;
data{3} = importdata('300Hz.txt') ;
data{4} = importdata('400Hz.txt') ;
data{5} = importdata('500Hz.txt') ;
data{6} = importdata('600Hz.txt') ;
data{7} = importdata('700Hz.txt') ;
data{8} = importdata('800Hz.txt') ;
data{9} = importdata('10kHz.txt') ;


% inddeler al data

for i = 1:9
    tid{i} = data{i}.data(:,1) * 10e-4 ;
    V_ind{i} = data{i}.data(:,2) ;
    V_ud{i} = data{i}.data(:,3) ;
end

%% 

% plot af r? data

for i = 1:1
    figure(i)
    plot(tid{i},V_ind{i})
    hold on
    plot(tid{i},V_ud{i})
    hold off
end

%%

amplitude_start = 0.5 ;
omega_start = 400*pi ;
phi_start = 1 ;

loop_start = 2 ;
loop_stop = 2 ;

for i = loop_start:loop_stop

    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    f = fittype('amplitude * sin(omega * x + phi)','options',s);
    fitobj = fit(tid{i},V_ind{i},f);
    figure(i)
    plot(tid{1},V_ind{1})
    hold on 
    plot(fitobj)
    %hold off
    
    amplitude_ind(i) = fitobj.amplitude ;
    omega_ind(i) = fitobj.omega ;
    phi_ind(i) = fitobj.phi ;
    
end
%%
for i = loop_start:loop_stop

    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    f = fittype('amplitude * sin(omega * x + phi)','options',s);
    fitobj = fit(tid{i},V_ud{i},f);
    plot(tid{1},V_ud{1})
    hold on 
    plot(fitobj)
    hold off
    
    amplitude_ud(i) = fitobj.amplitude ;
    omega_ud(i) = fitobj.omega ;
    phi_ud(i) = fitobj.phi ;
    
end
%%

plot(omega_ind,phi_ud-phi_ind,'o')

%%

V_ind = [9.0 ] ;
V_ud = [ ] ;
phi = [ ] ;