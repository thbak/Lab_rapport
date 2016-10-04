clc
clear all
close all

R = 9947; C = 47e-9;

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

data{10} = importdata('U_C.txt');
data{11} = importdata('U_R.txt');

% inddeler al data

for i = 1:11
    tid{i} = data{i}.data(:,1) * 10e-4 ;
    V_ind{i} = data{i}.data(:,2) ;
    V_ud{i} = data{i}.data(:,3) ;
end


%%

plot(tid{10},V_ind{10})
hold on
plot(tid{10},V_ud{10})
plot(tid{11},V_ud{11})

xlabel('t [s]','Fontsize',20)
ylabel('V [V]','Fontsize',20)
title('V_{in},V_C & V_R','Fontsize',20)
legend('V_{in}','V_C','V_R')


%% V_ind U_R U_C
interval = 31252:39244;

plot(tid{10}(interval),V_ind{10}(interval))
hold on
plot(tid{11}(interval),V_ud{11}(interval))
eps_start = 4;
RC_start = R*C; 

t = tid{10}(interval(end)) - tid{10}(interval(1));
A = mean(V_ind{10}(interval))*(1-exp( -1/(R*C)*t));

X = linspace(0,1.4e-3,1e3);
Y = A * exp( -1/(R*C).*X );
plot(X,Y)
legend('V_{ind}','U_R','teoretical','Location','Northeast');

f = fit(tid{11}(interval),V_ud{11}(interval),'exp1');
hold on
plot(f,'black')
disp('CR estimat')
(f.b)^(-1)
disp('CR teoretisk')
C*R
%%

amplitude_start_ind = 1/2 * ones(1,10);
amplitude_start_ud = 1/2 * [0.43 0.59 0.76 0.8 0.85 0.95 0.94 1 1 1.1];
omega_start_ind = 2*pi* [108 201 300 339 400 500 599 699 901 999];
phi_start_ud = 10^(-3)*omega_start_ind.* [8.1 4.2 2.9 2.6 2.2 1.8 1.5 1.3 1.0 0.9];

loop_start = 1 ;
loop_stop =  9 ;

for i = loop_start:loop_stop

    amplitude_start = amplitude_start_ind(i);
    omega_start = omega_start_ind(i);
    phi_start = 0 ;
    
    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    f = fittype('amplitude * sin(omega * x + phi)','options',s);
    fitobj = fit(tid{i},V_ind{i},f);
    fig = figure(1);    set(fig,'Position',[300 150 800 500]);
    subplot(3,3,i)
    plot(tid{i},V_ind{i},'color',[0 0 0])
    xlabel('t')
    ylabel('V')
    hold on 
    plot(fitobj)
    xlabel('t')
    ylabel('V')
    legend('V_{ind}','fit')
    hold off
    
    amplitude_ind(i) = fitobj.amplitude ;
    omega_ind(i) = fitobj.omega ;
    phi_ind(i) = fitobj.phi ;
    
end
%%



for i = loop_start:loop_stop
    
    amplitude_start = amplitude_start_ud(i);
    omega_start = omega_start_ind(i);
    phi_start = phi_start_ud(i)-phi_ind(i) ;
    
    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    %name=strcat('amplitude*sin(',num2str(omega_start),'*x+phi)');
    f = fittype('amplitude*sin(omega*x+phi)','options',s);
    fitobj = fit(tid{i},V_ud{i},f);
    fig = figure(2);
    set(fig,'Position',[300 150 800 500]);
    subplot(3,3,i)
    plot(tid{i},V_ud{i})
    hold on 
    plot(fitobj)
    xlabel('t')
    ylabel('V')
    legend('V_{ud}','fit')
    hold off
    
    amplitude_ud(i) = fitobj.amplitude ;
    omega_ud(i) = fitobj.omega ;
    phi_ud(i) = fitobj.phi ;
    
end
%% tilpasning af amplitude/phi

for i=1:9
    if phi_ud(i) >= 2*pi
        phi_ud(i) = phi_ud(i) - 2*pi;
    end
    if amplitude_ud(i) < 0
        phi_ud(i) = phi_ud(i) - pi;
        amplitude_ud(i) = abs(amplitude_ud(i));
    end
end

%% data_plot phi(wRC)

figure()
err = [           (0.011592783149663301)
                  (0.018646153184758468)
                  (0.021112905544267944)
                  (0.0210017309704781)
                  (0.01979647730448882)
                  (0.01826925924224841)
                  (0.016746549862393764)
                  (0.015299851443882018)
                  (0.01297396554327646)
];
errorbar(omega_ind*R*C,(phi_ud-phi_ind),err,'o','color',[0 0 0]);
%plot(omega_ind*R*C,(phi_ud-phi_ind),'o','color',[0 0 0])
axis([0 3 0 2])

% teoretisk phi(wRC):

hold on
wRC = linspace(0,3,1000);
g = atan(1./(wRC));
plot(wRC,g,'color',[0 0 0]);
xlabel('\omega RC [1]','Fontsize',20)
ylabel('\phi [1]','Fontsize',20)
title('\phi as a function of \omegaRC','Fontsize',20)
axis([0 3 0 2])
box on



%% chi^2 for phi(wRc)
chi2 = sum( (phi_ud - atan(1./(R*C.*omega_ind) )).^2 ./ (2*pi*0.1)^2);


%% plot af V_ind & V_ud

figure()
plot(omega_ind,amplitude_ind)
hold on
%plot(omega_ind,amplitude_ud,'o')
X = linspace(2*pi*100,2*pi*1000,1000);
Y = mean(amplitude_ind)*R ./ (sqrt(R.^2+(1./(X.*C)).^2));
err = [                                    (0.005591815974901736)
                  (0.008159317166617122)
                  (0.008174650969791645)
                  (0.007207095079059843)
                  (0.006127340959962654)
                  (0.005238762565048534)
                  (0.004574195569680601)
                 (0.0040884493341492014)
                 (0.0035381649975501613)

];
errorbar(omega_ind,amplitude_ud,err,'o','color',[0 0 0])
plot(X,Y,'color',[0 0 0])
xlabel('\omega','Fontsize',20)
ylabel('Voltage [V]','Fontsize',20)
title('V_{in} & V_{out} as a function of \omega','Fontsize',20)
legend('V_{in}','V_{out}','theoretical','Location','southeast')

A_start = mean(amplitude_ind);
C_start = C;

% s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[A_start,C_start]);
%     f = fittype('A*9947 / (sqrt(9947^2+(1/(x*C))^2))','options',s);
%     fitobj = fit(omega_ind,amplitude_ud,f);
%% teoretisk g(w)

w = linspace(0,1e4,1e3);
g = R ./ (sqrt(R.^2+(1./(w.*C)).^2));
figure()
plot(w,g,'color',[0 0 0])
xlabel('\omega [s^{-1}]','FontSize',20)
ylabel('g [1]','FontSize',20)
title('g as a function of \omega','FontSize',20)
box on

%% LAVPASFILTER!!!!!!!!!!!!!!!!!

R = 9947; C = 47e-9;

% importerer al data
data{1} = importdata('l9Hz.txt') ;
data{2} = importdata('l171Hz.txt') ;
data{3} = importdata('l400Hz.txt') ;
data{4} = importdata('l1000Hz.txt') ;
data{5} = importdata('l1800Hz.txt') ;
data{6} = importdata('l300Hz.txt') ;
data{7} = importdata('l700Hz.txt') ;
data{8} = importdata('l900Hz.txt') ;
data{9} = importdata('l100Hz.txt') ;

for i = 1:9
    tid{i} = data{i}.data(:,1) * 10e-4 ;
    V_ind{i} = data{i}.data(:,2) ;
    V_ud{i} = data{i}.data(:,3) ;
end

%%

amplitude_start_ind = 1/2 * ones(1,10);
amplitude_start_ud = 1/2 * [0.992 0.919 0.6746 0.325 0.199 0.7698 0.452 0.3733 0.9788];
omega_start_ind = 2*pi* [9.046 171 401.8 1052 1837 304.8 703.1 898.8 103];
phi_start_ud = 10^(-3)*omega_start_ind.* [0 0.4846 0.3448 0.1961 0.1284 0.3881 0.2622 0.222 0.4665];

loop_start = 1 ;
loop_stop =  9 ;

for i = loop_start:loop_stop

    amplitude_start = amplitude_start_ind(i);
    omega_start = omega_start_ind(i);
    phi_start = 0 ;
    
    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    f = fittype('amplitude * sin(omega * x + phi)','options',s);
    fitobj = fit(tid{i},V_ind{i},f);
    fig = figure(1);    set(fig,'Position',[300 150 800 500]);
    subplot(3,3,i)
    plot(tid{i},V_ind{i},'color',[0 0 0])
    xlabel('t')
    ylabel('V')
    hold on 
    plot(fitobj)
    xlabel('t')
    ylabel('V')
    legend('V_{ind}','fit')
    hold off
    
    amplitude_ind(i) = fitobj.amplitude ;
    omega_ind(i) = fitobj.omega ;
    phi_ind(i) = fitobj.phi ;
    
end

%%

for i = loop_start:loop_stop
    
    amplitude_start = amplitude_start_ud(i);
    omega_start = omega_start_ind(i);
    phi_start = phi_start_ud(i)-phi_ind(i) ;
    
    s = fitoptions('Method','NonLinearLeastSquares','Startpoint',[amplitude_start,omega_start,phi_start]);
    %name=strcat('amplitude*sin(',num2str(omega_start),'*x+phi)');
    f = fittype('amplitude*sin(omega*x+phi)','options',s);
    fitobj = fit(tid{i},V_ud{i},f);
    fig = figure(2);
    set(fig,'Position',[300 150 800 500]);
    subplot(3,3,i)
    plot(tid{i},V_ud{i})
    hold on 
    plot(fitobj)
    xlabel('t')
    ylabel('V')
    legend('V_{ud}','fit')
    hold off
    
    amplitude_ud(i) = fitobj.amplitude ;
    omega_ud(i) = fitobj.omega ;
    phi_ud(i) = fitobj.phi ;
    
end

%%