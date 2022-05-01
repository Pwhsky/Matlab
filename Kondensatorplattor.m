%% Numerisk uppgift – 2D fulla problemet

Nx = 50; % antalet funktionspunkter i x-led

Ny = 50; % antalet funktionspunkter i y-led

Ni = 80; % antalet iterationer

mpx = ceil(Nx/2); % mittpunkten i x-led, avrundat uppat

mpy = ceil(Ny/2); % mittpunkten i y-led, avrundat uppat



a = 10; % avstand fran sidor till plattor i x-led [cm]

b = 10; % avstand fran sidor till plattor i y-led [cm]

c = 1; % plattornas bredd [cm]

d = 2; % avstand mellan plattorna [cm]

e = 20; % plattornas langd [cm]



hx = (2*a+2*c+d)/Nx; % avstand mellan punkter i x-led, steglangd

hy = (2*b+e)/Ny; % avstand mellan punkter, steglangd

 

ap = floor(a/hx);

bp = floor(b/hy);

cp = floor(c/hx);

dp = floor(d/hx); 

ep = floor((e/hy)/2);



V1 = 5; % potentialen for platta 1

V2 = -5; % potentialen for platta 2

Vx = 0; % randpotentialen hogst upp & langst ner

Vy = 0; % randpotentialen till hoger & vanster



V = zeros(Nx,Ny); % potentialmatrisen



% Bestammer cellens randvillkor

V(1,:) = Vy;

V(Nx,:) = Vy;

V(:,1) = Vx;

V(:,Ny) = Vx;



% Bestammer potentialen i cellens horn

V(1,1) = 0.5*(V(1,2)+V(2,1));

V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));

V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));

V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));



% Plattornas position

posx = ceil(dp/2); % plattornas avstand fran mittpunkt i x-led

pp11 = mpx+posx;

pp12 = mpx+posx+cp;

pp21 = mpx-posx-cp;

pp22 = mpx-posx;





for z = 1:Ni 

        

        for i=2:Nx-1

        for j=2:Ny-1      

            

        %De tva forsta raderna gor så att plattornas potential ej andras

                V(pp11:pp12,mpy-ep:mpy+ep) = V1;

                V(pp21:pp22,mpy-ep:mpy+ep) = V2;

                

                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));

        end

        end

        

end



V = V';



[Ex,Ey]=gradient(V,hx,hy);

Ex = -Ex;

Ey = -Ey;



E = sqrt(Ex.^2+Ey.^2);  



x = (1:Nx)-mpx;

y = (1:Ny)-mpy;



K = [V(ceil(Nx/2)-ceil(2/hx),ceil(Ny/2)-ceil(2/hy)); 

    V(ceil(Nx/2)-ceil(2/hx),ceil(Ny/2)+ceil(2/hy));

    V(ceil(Nx/2)+ceil(2/hx),ceil(Ny/2)+ceil(2/hy));

    V(ceil(Nx/2)+ceil(2/hx),ceil(Ny/2)-ceil(2/hy))]; %konvergerande värden

%% Plotta potentialen

figure(1)

contour_range_V = -101:0.5:101;

contour(x,y,V,contour_range_V,'linewidth',0.5);

axis([min(x) max(x) min(y) max(y)]);

colorbar('location','eastoutside','fontsize',14);

xlabel('x-axel [cm]','fontsize',14);

ylabel('y-axel [cm]','fontsize',14);

title('Potentialen V(x,y) [V]','fontsize',14);

h1=gca;

set(h1,'fontsize',14);



% Plotta elektriska fältet

figure(2)

contour_range_E = -20:0.05:20;

contour(x,y,E,contour_range_E,'linewidth',0.5);

axis([min(x) max(x) min(y) max(y)]);

colorbar('location','eastoutside','fontsize',14);

xlabel('x-axel [cm]','fontsize',14);

ylabel('y-axel [cm]','fontsize',14);

title('Elektriska fältet E(x,y) [V/cm]','fontsize',14);

h2=gca;

set(h2,'fontsize',14);



% Plotta elektriska fältlinjerna

figure(3)

contour(x,y,E,'linewidth',1);

hold on, quiver(x,y,Ex,Ey,1)

title('Fältlinjer till E(x,y) [V/m]','fontsize',14);

axis([min(x) max(x) min(y) max(y)]);

colorbar('location','eastoutside','fontsize',14);

xlabel('x-axel [cm]','fontsize',14);

ylabel('y-axel [cm]','fontsize',14);

h3=gca;

set(h3,'fontsize',14);