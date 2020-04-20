% Heat conduction
 
clear; close all; clc;
n = 15;                %no. of divisions
T = ones(n,1); 
%Fourier's law

Temp = T;
T(1) = 10;          %end surface temperatures
T(n) = 1; 
x = linspace(0,1,n);    %length in m
dx = x(2)-x(1);
dt = dx.^2/5;
kappa = 1;       %thermal diffusivity
 
error = 1;
min = .001;
k = 0;
while error > min,
   Temp = T;
   k = k+1;
   for i = 2:n-1
    T(i) = dt*kappa*(Temp(i+1)-2*Temp(i)+Temp(i-1))/dx^2+Temp(i);
   end
   error = max(abs(T-Temp));
   if mod(k,13)==0, 
       y(k,:) = T; 
       plot(x,y)
        xlabel('x'),ylabel('Temperature'),
        title(['Heat Conduction']),
        pause(.2),
   end  
end

 
