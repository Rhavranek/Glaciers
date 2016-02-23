%Rachel Havranek
%February 2016, Geomodeling 2016

%% Initialize

%horizontal array
dx = 100;
xmax = 10000;
xmin=0;
x = 0:dx:xmax;

%slope array
zmax = 1000;
S = 0.1; %slope of the  bottom
zb = zmax - (S*x);

%initial ice:
z=zmax*ones(size(x));
H=1;
H=H*ones(size(x));

% time array
tmax=1000;
dt=.005; % it has to be really really small for the code to run
t=0:dt:tmax; %creates an array of time steps (time)
imax=length(t);
period=200;
nplots = 100; %so our code runs faster - we'll plot one hundred times
tplot = tmax/nplots; %so our code runs faster - divides total time up into 100 time steps

%Variables:
gamma=0.01; %per year
uslide=.05; %m/yr
rho=917; %kg/m^3
g=9.8; %m/s^2
A=2.1e-16; %pa
ELA=800;

%erosion:
e=0.02; %erosion constant 
erosion=e*uslide; %m/yr

%% RUN
for i=1:imax

%erosion from the ice
subice=find(z>zb); 
subair=find(z<=zb);
ei(subice)=erosion; %erosion occurs under the ice
ei(subair)=0; %but not under the air

%fluctuating ELA 
ELA=(400*sin((t(i)*2*pi)/period))+300; %sinusoidal ELA  

%glaciers loop
zb=zb-ei; %height of the bedrock after erosion
dzdx=diff(z)/dx; %change in height of the total (bedrock and height of ice)
b=gamma*(z-ELA); %ablation
Hedge=H(1:end-1) + (diff(H)/2);
Q=(uslide.*Hedge)+(A*((rho*g*abs(dzdx)).^3)).*((Hedge.^5)/5); %flux law
Q=[0 Q 0]; %no flux at the edges
dQdx=diff(Q)/dx; %change in flux
dHdt=b-dQdx; %change in height of ice
H=H +(dHdt.*dt); %total height of ice
H=max(0,H); %ice only exists above the surface
z=(zb+H); %total height of ice+bedrock


if (rem(t(i),tplot)==0)
figure(1)
plot(x,z,'k','linewidth',1.5)
     hold on
     plot(x,zb,'r','linewidth',1.5)
     hold off
      axis ([0,10000,0, 1500])
      xlabel('Distance (m)','fontname','arial','fontsize',21)
      ylabel('Height (m)','fontname','arial','fontsize',21)
     pause(0.2)
    
     
end
end


