# GitHub-Test

clear all;
clc
% One Dimensional Unsteady Heat Transfer Simulation Using Explicit Finite Difference Method
% Matlab using to solve numerically a heat transfer
% Leila Karimi, New Mexico State University, November 2011.
 
% L                       -m        Length of 1D space
% n                       -""       Number of nodes
% del_x                   -m        Spacing of nodes
% t0                     -Seconds  Starting time
% tmax                   -Seconds  Ending time
% k                       -W/(m K)  Thermal conductivity - conduction ciefficient
% alpha                   -m^2/s    Diffusivity
% h                       -W/m^2 K  Convection coefficient
% del_t                   -hour  Change in time each iteration
% T1                     -hour  Initial time
%h_out       average combined heat transfer coefficient for heat loss from
%the Trombe wall to the ambient is determined to be h_out=0.7 Btu/(h*ft^2F)
%h_in        heat transfer coefficient at the interior surface of the
%Trombe wall is h_in=1.8 Btu/(h*ft^2)
%%%% Trombe Wall Properties %
 
 
L= input(' Please input the thickness of the wall, (L=1) = ');
del_x= input(' Please input the uniform nodal spacing, (del_x=0.2) = ');
n=L/del_x+1;
%n=6;
T1= input(' Please input the internal surface of the wall, (T1(f)=70) = ');
%T1=70
T(n)= input(' Please input the external surface of the wall, (T(n)(f)=30) = ');
%T(n)=30;
 
%del_x=0.2;
k= input(' Please input conductivity coefficient,(K=0.4),= ');
%k=0.4;
t0= input(' Please input the minimum time (h), for this problem t0(h),(t0=7) = ');
%t0=7;
tmax = input(' Please input the maximum time (h), for two days=t0+48, (55) = ');
%tmax=55;
alpha_1 = input('Please input the thermal diffusivity, a(ft^2/s), (?= 4.78 x 10-6 ft2/s) = ');
alpha=alpha_1*3600;
%L=xmax-xmin;
 
%del_x= L/(n-1);
A=250;
%%%%%%%%   del_t should be less than 0.5*del_x^2/alpha   %%%%%%%%
del_t = 0.2*del_x^2/alpha ;
T_out=30;
T_in=70;
h_in= input(' Please input the heat transfer coefficient of internal fluid, h_in (1.8)  = ');
%h_in=1.8;
h_out=input(' Please input the heat transfer coefficient of external fluid, h_out (0.7)  = ');
%h_out=0.7;
Fo=alpha*del_t/(del_x^2);
B_in=h_in*del_x/k;
B_out=h_out*del_x/k;
 
ii=1
t_total=t0:del_t:tmax;
T_out=zeros(ii,6)
q=zeros(ii,6)
for ii=1:length(t_total)
    t = t_total(ii);
if t0<t & t<(t0+3)
     T_out(ii,6)=33 ; q(ii,6)=0.77*114;
    else if (t0+3)<=t & t<(t0+6)
    T_out(ii,6)=43 ; q(ii,6) =0.77*242;
    elseif (t0+6)<=t & t<(t0+9)
    T_out(ii,6)=45 ; q(ii,6) =0.77*178;
    elseif (t0+9)<=t & t<(t0+12)
    T_out(ii,6)=37 ; q(ii,6) =0;
    elseif (t0+12)<=t & t<(t0+15)
    T_out(ii,6)=32 ; q(ii,6) =0;
    elseif (t0+15)<=t & t<(t0+18)
    T_out(ii,6)=27 ; q(ii,6) =0;
    elseif (t0+18)<=t & t<(t0+21)
    T_out(ii,6)=26 ; q(ii,6) =0;
    elseif (t0+21)<=t & t<(t0+24)
    T_out(ii,6)=25 ; q(ii,6) =0;
    elseif (t0+24)<=t & t<(t0+27)
    T_out(ii,6)=33 ; q(ii,6)=0.77*114;
    else if (t0+27)<=t & t<(t0+30)
    T_out(ii,6)=43 ; q(ii,6) =0.77*242;
    elseif (t0+30)<=t & t<(t0+33)
    T_out(ii,6)=45 ; q(ii,6) =0.77*178;
    elseif (t0+33)<=t & t<(t0+36)
    T_out(ii,6)=37 ; q(ii,6) =0;
    elseif (t0+36)<=t & t<(t0+39)
    T_out(ii,6)=32 ; q(ii,6) =0;
    elseif (t0+39)<=t & t<(t0+42)
    T_out(ii,6)=27 ; q(ii,6) =0;
    elseif (t0+42)<=t & t<(t0+45)
    T_out(ii,6)=26 ; q(ii,6) =0;
    elseif (t0+45)<=t & t<(t0+48)
    T_out(ii,6)=25 ; q(ii,6) =0;
        end
    end
end
end
T_out;
q;
for ii=0
    for i=1:n
    T(1,i)=T1-((T1-T(n))/(n-1))*(i-1);
    end
 
end
T
 
t=t_total;
for ii=1:length(t)
    
    for i = 1:n
  
        if i == 1 
            
            T(ii+1,i) = (1-2*B_in*Fo-2*Fo)*T(ii,1)+2*Fo*(T(ii,i+1)+B_in*T_in);
            
            
        
        elseif i == n
            
       T(ii+1,i) = 2*Fo*(B_out*T_out(ii,6)+T(ii,n-1)+q(ii,6)*del_x/k)+(1-2*Fo-2*B_out*Fo)*T(ii,i);
      
        else
 
       T(ii+1,i) = Fo*(T(ii,i-1)+T(ii,i+1))+(1-2*Fo)*T(ii,i);
               
      
        end
 
end
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Calculating the node temperatures after 12, 24, 36 and 48 hours
 
    
T
for i=1:6
Temp12(1,i)=T(length(t)/4,i);
Temp24(1,i)=T(length(t)/2,i);
Temp36(1,i)=T(78,i);
Temp48(1,i)=T(length(t),i);
end
 
Temp12
Temp24
Temp36
Temp48

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the node temperatures after 12, 24, 36 and 48 hours
 
 
Temp=[Temp12;Temp24;Temp36;Temp48];
 
i=[1:1:6];
 
plot(i,Temp)
xlabel('node number ')
ylabel('Temperature (.F)')
title('Temperature changes for nodes after 12, 24, 36 and 48 hours')
legend('Temp12','Temp24','Temp36','Temp48',1)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating the net heat transfer to the home
%heat1 is for the first and heat2 is for the second day

heat1=zeros(1,length(t/2));
M=zeros(1,length(t/2));
for ii=1:length(t)/2
   heat1(1,ii)=ii;
M(1,ii)=h_in*A*(T(ii,1)-T_in);
          
end
heat1=trapz(heat1,M)
 
 
heat2=0;
for ii=length(t)/2:length(t)
     heat2(1,ii)=ii;
M(1,ii)=h_in*A*(T(ii,1)-T_in);
          
end
heat2=trapz(heat2,M)

%%%%5
% Plotting the  temperatures distribution in the form of surface for whole
% the two days continously along the wall
 
figure(1)
 
 hold on
 
figure(2)
colormap hot  
plot(T)
surf(T); figure (gcf)
xlabel('i (Node Number (Position through Wall)')
ylabel('ii (Time)')
zlabel('T (.F)')
title('Temperature profile changes within 48 hr')
 
shading interp
colorbar
