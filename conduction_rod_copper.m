clear;clc;
%-------------------------INPUT-----------------------------%
% number of divisions made
    N = 150;
% number of time steps
    time_steps = 50000;
% length of rod
    L = 0.250;
% radius of rod
    rad_1 = 0.01;
    rad = rad_1 .* ones(1,N);    
% thermal conductivity (in W/(m.degC)
    k_1 = 400;
    k = k_1 .* ones(1,N);
% specific heat capacity (in J/(kg.degC)
    c_1 = 380;  
    c = c_1 .* ones(1,N);
% density  (in kg/m^3)
    rho_1 = 8900;
    rho = rho_1 .* ones(1,N);
% initial temperature of rod
    T = zeros(time_steps,N);
% intial energy flux density
    J = zeros(time_steps,N+1);

%-----------------BOUNDARY CONDITIONS-----------------------%  

% End surface temperatures
    T(:,1) = 100;
    T(:,end) = 0;
    

%---------------------CALCULATIONS-----------------------------%    

dx = L / N                 % length of elements
dt = (10/3)*dx              % increment in time
x = dx/2 + dx .*(0:N-1);
xi = 0:dx:L;
t = 0 : (time_steps-1);
t = dt .* t;
K_1 = (k ./ dx);
K_2 = dt ./ (rho .* c .* dx);
area = pi * rad(1)^2;       % Area of cross section

for i = 1 : time_steps-1
for j = 2 : N
    J(i+1,j) = K_1(j) * ( T(i,j-1) - T(i,j) );
    J(i+1,1) = J(i+1,2);
    J(i+1,end) = J(i+1,end-1);
end  

for j = 2 : N-1
    T(i+1,j) = T(i,j) + K_2(j) * ( J(i+1,j) - J(i+1,j+1) ) ;
end
end

% energy flux dQ/dt = JA
    dQdt = J(end,end) .* area
% temperature gradient dT/dx
    dTdx = (T(end,end) - T(end,1)) / (x(end) - x(1))
% max time
t_Max = t(end)

%-----------------------------Plots-----------------------------%

                        figure(1)
 %------------------Energy flux density vs position-------------------%
 index = round(0.1 * time_steps) .* [1/(0.1*time_steps) 1 2 5 10];
 x_plot = xi;
 y_plot = J(index,:);
   plot(x_plot,y_plot);
   xlabel('position  x  in [m]');
   ylabel('flux density  J in [w/m^2]');
   h_title = title('time elapsed in [seconds]')
   time_1 = '0';
   time_2 = num2str(t(index(2)));
   time_3 = num2str(t(index(3)));
   time_4 = num2str(t(index(4)));
   time_5 = num2str(t(index(5)));
   h_legend = legend(time_1,time_2,time_3,time_4, time_5);

figure(2)
   x_plot = x;
   y_plot = T(index,:);
   plot(x_plot,y_plot);
   xlabel('position  x in [m]');
   ylabel('Temperature T in [degC]');
   
                        figure(3)
%----------------- energy flux density J vs time plot -----------------%
index = round((0.25 * N) .* [1/(0.25*N) 1 2 3 4]);

   x_plot = t; y_plot = J(:,index);
   plot(x_plot,y_plot);
   xlabel('time  t  [s]');
   ylabel('flux density  J  [w/m^2]');
   h_title = title('x position along rod in meters');
   
   position_1 = '0';
   position_2 = num2str(x(index(2)));
   position_3 = num2str(x(index(3)));
   position_4 = num2str(x(index(4)));
   position_5 = num2str(x(index(5)));
   h_legend = legend(position_1,position_2,position_3,position_4, position_5);

                        figure(4)
%-----------------temperature T vs time t-------------------%
   x_plot = t; y_plot = T(:,index);
   plot(x_plot,y_plot);
   xt = 'time  t  [s]';
   yt = 'Temperature T  [degC]';
   xlabel('time  t in [s]');
   ylabel('Temperature T in [degC]');