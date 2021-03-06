%Rachel Havranek
%February 2016, Geomodeling 2016

%% Initialize
clear all
figure(1)
clf

%horizontal array
dx = 200;
xmax = 30000;
xmin=0;
x = 0:dx:xmax;
xstar=2000;

%slope array
zmax = 3000;
zmax2=800;
zmin=1600;
%zb=zmin+((zmax-zmin)*exp(-x/xstar));
S = 0.1; %slope of the  bottom - similar to Macgregor
z1 = zmax - (S*x);
z2=zmax2*exp(-x./xstar);
zb=z1+z2;

%initial ice:
z=zmax*ones(size(x));
H=1;
H=H*ones(size(x));

% time array
tmax=2000;
dt=.002; % it has to be really really small for the code to run
t=0:dt:tmax; %creates an array of time steps (time)
imax=length(t);
period=50;
nplots = 100; %so our code runs faster - we'll plot one hundred times
tplot = tmax/nplots; %so our code runs faster - divides total time up into 100 time steps

%Variables:
gamma=0.01; %per year
uslide=.05; %m/yr
rho=917; %kg/m^3
g=9.8; %m/s^2
A=2.1e-16; %pa
% ELA=800;
ELAplot=ones(size(x));

%erosion:
e=.03; %erosion constant - this sets my erosion at 1.5 mm/yr this fits w/ Macgregor
erosion=e*uslide; %m/yr 
nframe = 0;

%% RUN
for i=1:imax

%erosion from the ice
subice=find(z>zb); 
subair=find(z<=zb);
ei(subice)=erosion; %erosion occurs under the ice
ei(subair)=0; %but not under the air

%fluctuating ELA 
ELA=(400*sin((t(i)*2*pi)/period))+1900; %sinusoidal ELA  

%glaciers loop
zb=zb-ei; %height of the bedrock after erosion
dzdx=diff(z)/dx; %change in height of the total (bedrock and height of ice)
b=gamma*(z-ELA); %ablation
Hedge=H(1:end-1) + (diff(H)/2);
Q=(uslide.*Hedge)+(A*((rho*g*abs(dzdx)).^3)).*((Hedge.^5)/5); %flux law
Q=[0 Q Q(end)]; %no flux at the top edge, but the flux at the end of the array defines the end flux
dQdx=diff(Q)/dx; %change in flux
dHdt=b-dQdx; %change in height of ice
H=H +(dHdt.*dt); %total height of ice
H=max(0,H); %ice only exists above the surface
z=(zb+H); %total height of ice+bedrock
agl(i)=(sum(H))*dx;


if (rem(t(i),tplot)==0)
    nframe = nframe+1
figure(1)
plot(x/1000,z,'k','linewidth',1.5)
     hold on
     plot(x/1000,zb,'r','linewidth',1.5)
     above=find(z<ELA);
     plot(x(above)/1000,ELA*ones(size(x(above))),'b--','linewidth', 2)
     ht=text(25,3000,['  ',num2str(t(i)), 'years '],'fontsize',18);
      axis ([0,xmax/1000,-1000, 4000])
      xlabel('Distance (km)','fontname','arial','fontsize',21)
      ylabel('Height (m)','fontname','arial','fontsize',21)
      drawnow
      hold off
      
     pause(0.5)
         
end
end

  figure (2)
     plot (t,agl/1000000)
     hold off
     xlabel('Time','fontname','arial','fontsize',21)
     ylabel('Volume (km^3)','fontname','arial','fontsize',21)
     
  figure (3)
  plot (t, ELA)
  hold off
