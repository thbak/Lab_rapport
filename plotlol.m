data = importdata('speanding_over_resistor.txt') ;

tid = data.data(:,1) ;
V_ind = data.data(:,2) ;
V_r = data.data(:,3) ;

plot(tid,V_ind)
hold on
plot(tid,V_r)

omega = (47*10^(-9) * 10^4)^(-1) ;

%% 

data = importdata('omega_0.txt') ;

tid = data.data(:,1) ;
V_ind = data.data(:,2) ;
V_r = data.data(:,3) ;

plot(tid,V_ind)
hold on
plot(tid,V_r)

%% 

data = importdata('U_C.txt') ;

tid = data.data(:,1) ;
V_ind = data.data(:,2) ;
V_r = data.data(:,3) ;

plot(tid,V_ind)
hold on
plot(tid,V_r)

%% 

data = importdata('U_R.txt') ;

tid = data.data(:,1) ;
V_ind = data.data(:,2) ;
V_r = data.data(:,3) ;

plot(tid,V_ind)
hold on
plot(tid,V_r)