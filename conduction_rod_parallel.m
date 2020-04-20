clear;clc;
%-----------------------------INPUTS--------------------------------

% number of divisions made
    N = 150;
    n(1) = 70; n(2) = 20 + n(1);
% number of time steps
    time_steps = 70000;
% length of rods (in m)
    L = 0.25;
% radius of rod (in m)
    r1 = .01;
    rad = r1 .* ones(1,N);    
% thermal conductivity of rods (in W/(m.degC)
    k_1 = 400;
    k_2 = 350;
    k_3 = 325;
    k = k_1 .* ones(1,N);
    k(n(1)+1:end) = k_2;
    k(n(2)+1:end) = k_3;

% specific heat capacity of rods (in J/(kg.degC)
    c_1 = 380;  
    c_2 = 350;
    c_3 = 320;
    c = c_1 .* ones(1,N);
    c(n(1)+1:end) = c_2;
    c(n(2)+1:end) = c_3;

% density of rods (in kg/m^3)
    rho1 = 8900;
    rho2 = 8000;
    rho3 = 5000;
    rho = rho1 .* ones(1,N);
    rho(n(1)+1:end) = rho2;
    rho(n(2)+1:end) = rho3;

    
% initial temperature of rods  (in degC)
    T = zeros(time_steps,N);

% intial energy flux density
    J = zeros(time_steps,N+1);

%------------------ BOUNDARY CONDITIONS -------------------%
    T(:,1) = 100;
    T(:,end) = 0;
    T(:,1:end) = 100;
    T(:,end/2:end) = 0;
%------------------------ CALCULATIONS --------------------------%
dx = L / N;                  % width of elements
dt = (10/3)*dx                 %increment in time
x = dx/2 + dx .*(0:N-1);
xi = 0:dx:L;
t = 0 : (time_steps-1);
t = dt .* t;
K_1 = (k ./ dx);                            
K_2 = dt ./ (rho .* c .* dx);                
area = pi * rad(1)^2;               % cross-sectional area
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
    
%-----------------------------plots--------------------------------%
                            figure(1);
%-------------- energy flux density J vs position -----------------%
index = round(0.1 * time_steps) .* [1/(0.1*time_steps) 1 2 5 10];

   xplot = xi; 
   yplot = J(index,:);
   plot(xplot,yplot);
   xt = 'position  x  [m]';
   yt = 'flux density  J  [w/m^2]';
   xlabel(xt);
   ylabel(yt);
   h_title = title('Time in [s]');
   time_1 = '0';
   time_2 = num2str(t(index(2)));
   time_3 = num2str(t(index(3)));
   time_4 = num2str(t(index(4)));
   time_5 = num2str(t(index(5)));
   h_legend = legend(time_1,time_2,time_3,time_4, time_5);
   
                            figure(2)

% -----------------temperature T vs position ----------------------%
  index = round(0.1 * time_steps) .* [1/(0.1*time_steps) 1 2 5 10];
  xplot = x; yplot = T(index,:);
   xplot = x; yplot = T(index,:);
   plot(xplot,yplot);
   xt = 'position  x in [m]';
   yt = 'Temperature T in [degC]';
   xlabel(xt);
   ylabel(yt);