% 09:15 14/4-2018
% Sebastian Bertolino MATLAB R2014b
% FYP104 V�gr�relsel�ra (moment Programmering med matlab - 1.0 hp)
% Deadline Rapport inl�mning (f�r student-feedback): Onsdag 2:e Maj klockan 23:55
% Feedback -inl�mning: Onsdag 9:e Maj klockan 23:55
% Rapporten inl�mnas senast: Onsdag 16:e maj klockan 23:55
% Koden ska med i rapporten!

% Uppgift 1: Interferens fr�n ett antal spalter
N = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]; % Antal spalter
l = 10; % L�ngd fr�n spalter till sk�rmen
b = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avst�nd mellan spalterna
r = [-1,1]; % Utv�rdera |S(x)| i detta intervall i cm

% l >>  d, t.ex l = 10cm

% Plotta amplituden | S(x)| av den totala v�gen p� sk�rmen som funktion av avst�ndet fr�n
% centrumstr�len f�r n�gra olika v�rden p� N = 1, 2, 3, ... Beskriv hur interferensm�nstret
% �ndras n�r N �ndras. F�rs�k f�rklara varf�r. Ge om m�jligt ett generellt uttryck f�r detta?

% Plotta |S(x)| f�r n�gra olika v�gl�ngder (allt annat konstant), t.ex. r�tt ljus 700 nm, gr�nt 
% ljus 500 nm och bl�tt ljus 450 nm. Beskriv �terigen vad som h�nder n�r lambda �ndras och varf�r.
% �ndra p� avst�ndet mellan spalterna. Hur f�r�ndras |S(x)| n�r b minskar?
% Varf�r sker det?
% �ndra p� avst�ndet mellan spalterna. Hur f�r�ndras jS(x)j n�r b minskar? Varf�r sker det?
lambda = [700*10^(-5), 500*10^(-5), 450*10^(-5)]; % v�gl�ngden i cm 700: r�tt, 500: gr�nt, 450: bl�tt
s_s = r(1); % Star p� sk�rmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d = abs((r(1)-r(2))/res); % Distans mellan punkterna
S = zeros(size(s_s:s_d:s_e,2), size(lambda,2), size(b,2), size(N,2)); % Amplituden i alla punkter p� sk�rmen


% Vi ska nu r�kna ut amplituden, detta g�r vi med genom att r�kna alla
% bidrag fr�n spalterna i varje punkt.
% S(alpha,beta,gama) alpha: x p� sk�rmen, beta: lambda, gama: spaltavst�nd
n_count = 1;
for m = 1:length(lambda) % �ver alla lambda
    for k = 1:length(b) % �ver alla b, distansen mellan punkterna
        for t = 1:size(S,1) % �ver alla sk�rm punkter
            n_count=1;
            for N_temp = N
                for j = 1:N_temp % �ver alla spalter % Fel verkar komma fr�n -j*b(k) f�r den r�r sig med b, �kar med N ocks�
                    %S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
                    %S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-j*b(k)*(N_temp-1)/2)^2 ) /lambda(m)) );
                    S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-b(k)*(j-1-((N_temp-1)/2)))^2 ) /lambda(m)) );
                    %j-1 f�r att f� nollpunkten r�tt, f�r j startar som
                    %j=1.
                end
                % Vi m�ste sprida spalteran j�mt runt nollan, tror nu att
                % alla r�knar fr�n l�ngst ner. Eller? varfr� ser d�
                % spalterna med h�gt b s� bra ut? det borde ju vara tv�rtom
                % isf!
                n_count = n_count + 1;
            end
        end
    end
end
N_count = 1;
for N_tempplot = N
    % Pr�vade -100:0.1:100lechal@student.
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i r�knar alla i f�rsta
    % raden sen n�sta rad tills den har kommit till v�rdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b)
            subplot(length(lambda), length(b), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(s_s:s_d:s_e,2)-1)/2):((size(s_s:s_d:s_e,2)-1)/2), abs(S(:,i,j,N_count)))
            % n�r jag k�r den l�nga s� f�ljer bara halva med
            n=n+1;
        end
    end
    % Source
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    
    % Ger text till titlar och axlarna
    p=1;
    for i = 1:length(lambda)
        for j = 1:length(b)
            ax = findobj(f,'Type','Axes');
            title(ax(p),['b = ', num2str(b(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan s� va [] ist�llet {}, d� br�t den varje g�ng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d), ' [cm]' ])
            xlim(ax(p),[-(size(S,1)-1)/2-1 (size(S,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['Antal spalter: ' num2str(N(N_count)) ', l�ngd fr�n spalt till sk�rm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end

% Om man ska l�sa den generellt f�r varf�r den varierar med N, s� dela upp
% varje spalt f�r sig. Sen s� ser du hur m�nga spalt steg det tar f�r att f�
% motsatt i visa punkter. Borde finnas ett generellt uttryck f�r det.

%Solved:
% hur ska jag f� fram 0.001 fr�n -1 till 1 och 2000 punkter.
% abs(1-(-1)/2000 = 0.001, i algebra form: (r(1)-r(2))/

%%

% Uppgift 2: Interferens fr�n enkelspalt
% [-a/2, a/2],  S(x) = (1/a) Int_{-a/2}^{a/2} e^{i2Pi d(x,y) / ?} dy
% l�mpligt v�rde p� a ~ 10 \lambda, kan ocks� beh�vas att kolla p� x = [-3,3]

% Plotta intensiteten, d.v.s. kvadraten p� amplituden I = |S(x)|^2, p� interferensm�nstret
% fr�n en enkelspalt. Hur �ndras diffraktionsm�nstret n�r spaltbredden a eller v�gl�ngden
% ? �ndras?a(1);


%%% TEMP
if exist('lambda','var') == 0
    lambda = [700 *10^(-5), 500*10^(-5), 450*10^(-5)];
end

% Kollar om a existerar
if exist('a','var') ==1 
    a = 10*lambda(1);
end

if exist('r','var') == 0
    r = [-3 3];
end
if exist('l','var') == 0
    l = 10;
end

%%% TEMP

lambda_tmp = 100*lambda(1);
res = 200;
s_d_2 = abs((r(1)-r(2))/res);



%x = 0; % linspace(-3,3); % singel spalt placerad precis i nollpunkten
d = @(x,y) sqrt(l^2+(x-y).^2);

fun = @(x, y, a) (1/a)*exp(1i*2*pi*d(x,y)/lambda_tmp);  % fun(-a/2,-3,a)
% Varf�r f�r jag s� stroa v�rden? a �r litet men �nd�, bara exp delen ger
% 1.665e+38
 
% Hur skickar vi in x v�rdet fr�n integral funktionen?

% plot(-10:10, fun(-10:10, 2, 10*lambda_tmp, lambda_tmp)
% varf�r f�r jag s� stora v�rden?

%(1/10*lambda_tmp)*integral(@(x)exp(i*2*pi*d(x,y)/lambda_tmp),-3,3,'ArrayValued',true);
% integrera �ver alla sm� delar av spalten
fun2 = @(x)integral(@(y)fun(x, y, a), -a/2, a/2, 'ArrayValued', true);
% Fr�n uppgiften ovan: exp( (1i*2*pi*sqrt(l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
% 
close all %stänger alla figurer s� att mega plottningen nedanför funkar och vi inte döpper om massa irrelevanta saker.
% summera i alla punkter x
S_2 = zeros(size(r(1):s_d_2:r(2),2),1);
for i = 1:size(S_2)
    S_2(i)=fun2(i*s_d_2+r(1));
end
figure(98)
plot(r(1):s_d_2:r(2), abs(S_2).^2)
title('Intensitet fr�n enkel spalt')
xlabel('x [cm]')
ylabel('Amplitud')

% Nu ska vi kolla alla variationer av v�rden.
b_2 = a*[0.01 0.1 0.5 1 2 10 100]; % b �r variationen av a
% u spaltbrädd variabel, k �r lambda
fun = @(x, y, u, k) (1/b_2(u))*exp(1i*2*pi*d(x,y)/lambda(k));

fun2 = @(x, k, u)integral(@(y)fun(x, y, u, k), -b_2(u)/2, b_2(u)/2, 'ArrayValued', true);

S_2 = zeros( size(r(1):s_d_2:r(2), 2), size(lambda, 2), size(b_2, 2) );

for j = 1:size(S_2,1) % x
    for k = 1:length(lambda) % v�gl�ngden
        for u = 1:length(b_2) % spalt
            temp=fun2(j*s_d_2+r(1), k, u);
            S_2(j,k,u) = temp;
        end
    end
end

N_count = 1;
N_2=1;
for N_tempplot = N_2
    % Pr�vade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i r�knar alla i f�rsta
    % raden sen n�sta rad tills den har kommit till v�rdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            subplot(length(lambda), length(b_2), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_2:r(2),2)-1)/2):((size(r(1):s_d_2:r(2),2)-1)/2), abs(S_2(:,i,j)).^2)
            % n�r jag k�r den l�nga s� f�ljer bara halva med
            n=n+1;
        end
    end
    % Source
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    
    % Ger text till titlar och axlarna
    p=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            ax = findobj(f,'Type','Axes');
            title(ax(p),['b = ', num2str(b_2(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan s� va [] ist�llet {}, d� br�t den varje g�ng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_2), ' [cm]' ])
            xlim(ax(p),[-(size(S_2,1)-1)/2-1 (size(S_2,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['L�ngd fr�n spalt till sk�rm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count + 1;
    lambda = [700*10^(-5), 500*10^(-5), 450*10^(-5)]; % v�gl�ngden i cm 700: r�tt, 500: gr�nt, 450: bl�tt
s_s = r(1); % Star p� sk�rmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d = abs((r(1)-r(2))/res); % Distans mellan punkterna1;
end



% vad menar han med kontinuum gr�nsen?


% Tips: Definera funktionen d(x,y) = sqrt(l^2 + (x-y)^2 ) som en anonym
% funktion i matlab d�r x_j -> y �r en kontinuumgr�nsen av spalterna i
% f�reg�ende uppgift och x �r som f�rut. Integrera med matlabfunktionen
% integral. ( Man kan imtegrera em funktion av tv� variabler med avseende
% p� endast en p� f�ljande vis: fun2 =
% @(a)integral(@(b)fun(a,b),1,2,'ArrayValued',true);

% Source: https://se.mathworks.com/matlabcentral/answers/169362-check-if-variable-exists-in-workspace-to-plot-variable-else-generate-error
% Kollar om lambda existerar
if exist('lambda','var') == 0
    lambda = [700 *10^(-5), 500*10^(-5), 450*10^(-5)];
end

% Kollar om a existerar
if exist('a','var') ==1 
    a = 10*lambda(1);
else 
    a = 10*10*[700*10^(-5), 500*10^(-5), 450*10^(-5)];
end

l = 10;
r = [-3, 3];
a = 10*lambda(1);
lambda_tmp = 100*lambda(1);
res = 200;
s_d_2 = abs((r(1)-r(2))/res);

%x = 0; % linspace(-3,3); % singel spalt placerad precis i nollpunkten
d = @(x,y) sqrt(l^2+(x-y).^2);

fun = @(x, y, a) (1/a)*exp(1i*2*pi*d(x,y)/lambda_tmp);  % fun(-a/2,-3,a)
% Varf�r f�r jag s� stroa v�rden? a �r litet men �nd�, bara exp delen ger
% 1.665e+38
 
% Hur skickar vi in x v�rdet fr�n integral funktionen?

% plot(-10:10, fun(-10:10, 2, 10*lambda_tmp, lambda_tmp)
% varf�r f�r jag s� stora v�rden?

%(1/10*lambda_tmp)*integral(@(x)exp(i*2*pi*d(x,y)/lambda_tmp),-3,3,'ArrayValued',true);
% integrera �ver alla sm� delar av spalten
fun2 = @(x)integral(@(y)fun(x, y, a), -a/2, a/2, 'ArrayValued', true);
% Fr�n uppgiften ovan: exp( (1i*2*pi*sqrt(l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
% 
close all %stänger alla figurer så att mega plottningen nedanför funkar och vi inte döpper om massa irrelevanta saker.
% summera i alla punkter x
S_2 = zeros(size(r(1):s_d_2:r(2),2),1);
for i = 1:size(S_2)
    S_2(i)=fun2(i*s_d_2+r(1));
end
figure(98)
plot(r(1):s_d_2:r(2), abs(S_2).^2)

% Nu ska vi kolla alla variationer av värden.
b_2 = a*[0.01 0.1 0.5 1 2 10 100]; % b är variationen av a
% u spaltbrädd variabel, k är lambda
fun = @(x, y, u, k) (1/b_2(u))*exp(1i*2*pi*d(x,y)/lambda(k));

fun2 = @(x, k, u)integral(@(y)fun(x, y, u, k), -b_2(u)/2, b_2(u)/2, 'ArrayValued', true);

S_2 = zeros( size(r(1):s_d_2:r(2), 2), size(lambda, 2), size(b_2, 2) );

for j = 1:size(S_2,1) % x
    for k = 1:length(lambda) % våglängde
        for u = 1:length(b_2) % spalt
            temp=fun2(j*s_d_2+r(1), k, u);
            S_2(j,k,u) = temp;
        end
    end
end

N_count = 1;
N_2=1;
for N_tempplot = N_2
    % Pr�vade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i r�knar alla i f�rsta
    % raden sen n�sta rad tills den har kommit till v�rdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            subplot(length(lambda), length(b_2), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_2:r(2),2)-1)/2):((size(r(1):s_d_2:r(2),2)-1)/2), abs(S_2(:,i,j)).^2)
            % n�r jag k�r den l�nga s� f�ljer bara halva med
            n=n+1;
        end
    end
    % Source
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    
    % Ger text till titlar och axlarna
    p=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            ax = findobj(f,'Type','Axes');
            title(ax(p),['b = ', num2str(b_2(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan s� va [] ist�llet {}, d� br�t den varje g�ng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_2), ' [cm]' ])
            xlim(ax(p),[-(size(S_2,1)-1)/2-1 (size(S_2,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['L�ngd fr�n spalt till sk�rm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end



% vad menar han med kontinuum gr�nsen?

%%

% Uppgift 3: Interferens fr�n flera spalter av finit utbredning
% Studera nu hur interferensm�nstret �ndras n�r man tar h�nsyn till att varje spalt har en
% finit utbredningen enligt resultatet i uppgift 2. Plotta intensiteten I = |S(x)|^2 f�r n�gra
% olika spaltbredder. Beskriv hur interferensm�nstret fr�n flera enkelspalter ser ut i relation
% till de fr�n uppgift 1 och 2. V�lj i �vrigt alla parametrar s� som ovan

% variera a, 

if exist('N','var') == 0
    N = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]; % Antal spalter
end
if exist('l','var') == 0
    l = 10; % L�ngd fr�n spalter till sk�rmen
end
if exist('lambda','var') == 0
    lambda = [700*10^(-5), 500*10^(-5), 450*10^(-5)];
end
if exist('a','var') == 0
    a = lambda(1);
end
if exist('b','var') == 0
    b = a*[0.01 0.1 0.5 1 2 10 100];
    b_3 = b;
end
if exist('r','var') == 0
    r = [-1 1];
end
if exist('r','var') == 1
    s_s = r(1);
    s_e = r(2);
end
s_s = r(1); % Star p� sk�rmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d_3 = abs((r(1)-r(2))/res); % Distans mellan punkterna



% Nu ska vi kolla alla variationer av värden.
b_3 = a*[0.01 0.1 0.5 1 2 10 100]; % b �r variationen av a
% u spaltbrädd variabel, k är lambda
d_3 = @(x,y,k,j)sqrt( l^2 + (s_s+x-k*(j-1-((N_temp-1)/2)))^2 );
%d_3 = @(x,y,k,j)sqrt( l^2 + (s_s+t*s_d_3-k*(j-1-((N_temp-1)/2)))^2 );

fun = @(x, y, u, k) (1/(u))*exp(1i*2*pi*d_3(x,y,u,j)/lambda(k));

fun2 = @(x, k, u, lower)integral(@(y)fun(x, y, u, k), u, lower, 'ArrayValued', true);

S_3 = zeros( size(r(1):s_d_3:r(2), 2), size(lambda, 2), size(b_3, 2), size(N, 2) );

for j = 1:size(S_3,1) % x
    for k = 1:length(lambda) % våglängden
        for u = 1:length(b_3) % spaltbred
            n_count=1;
            for N_temp = N
                for g = 1:N_temp % �ver alla spalter % Fel verkar komma fr�n -j*b(k) f�r den r�r sig med b, �kar med N ocks�
                    
                    S_3(j,k,u,n_count) =fun2(j*s_d_3+r(1), k, b_3(u)*(g-(N_temp-1)/2)-a/2, b_3(u)*(g-(N_temp-1)/2)-a/2+a/2);
                    
                end
                % Vi m�ste sprida spalteran j�mt run = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avst�nd mellan spalternat nollan, tror nu att
                % alla r�knar fr�n l�ngst ner. Eller? varfr� ser d�
                % spalterna med h�gt b s� bra ut? det borde ju vara tv�rtom
                % isf!
                n_count = n_count + 1;
            end
            
        end
    end
end



N_count = 1;
N_3=[1 3 5 10 20 50 100];
for N_tempplot = N_3
    % Pr�vade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i r�knar alla i f�rsta
    % raden sen n�sta rad tills den har kommit till v�rdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_3)
            subplot(length(lambda), length(b_3), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_3:r(2),2)-1)/2):((size(r(1):s_d_3:r(2),2)-1)/2), abs(S_3(:,i,j)).^2)
            % n�r jag k�r = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avst�nd mellan spalterna den l�nga s� f�ljer bara halva med
            n=n+1;
        end
    end
    % Source
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    % https://se.mathworks.com/matlabcentral/answers/23204-how-to-give-labels-and-title-to-all-subplot-one-time
    
    % Ger text till titlar och axlarna
    p=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            ax = findobj(f,'Type','Axes');
            title(ax(p),['Antal spalter: ' num2str(N_3(N_count)) 'b = ', num2str(b_3(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan s� va [] ist�llet {}, d� br�t den varje g�ng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_3), ' [cm]' ])
            xlim(ax(p),[-(size(S_3,1)-1)/2-1 (size(S_3,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['L�ngd fr�n spalt till sk�rm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end




%%
% Bilda fouriertransformen av vektorn (med matlabfunktionen fft) och plotta dess absolutbelopp
% (d.v.s. den s� kallade spektralfunktionen) f�r n�gra olika antal spalter och olika avst�nd mellan spalterna.
% Finns det n�gra likheter med resultatet i Uppgift 1? F�rs�k beskriva vad fouriertransformen g�r,
% hur gittret och interferensm�nstret �r relaterat.
% Hur skulle man g�ra fft p� uppgift 2 och 3 f�r att f� motsvarande resultat f�r interferens
% fr�n en eller flera enkelspalter?

% Tips: Den numeriska algoritmen f�r att ber�kna Fourier-transformen �r mest effektiv om
% vektorn �r av l�ngden 2n med n ett heltal, n = 10 funkar utm�rkt.

% Tips: L�s hj�lpen f�r fftshift och fundera p� hur den han vara anv�ndbar.

spalt_N = [2 4 8 9 16 25 51 101]; % med olika space och spalt_N s� blir matrisen inte lika l�ng p� alla axlar.
space = [1 3 4 6 8];
space_slit = 1;

% borde g�ra det med cells.
spalt = zeros(max(space)*(max(spalt_N)+1)+max(space_slit)*max(spalt_N),size(spalt_N,2),size(space,2));

for j = 1:length(spalt_N)
    for k = 1:length(space)
        for i = 1:spalt_N(j)
            spalt((space(k)+1)*i,j,k) = 1;
        end
    end
end

figure(401)
plot(abs(fft(spalt(:,1,1))))

% suplot(x,y,i) x: antal rader, y: antal coloumner, i r�knar alla i f�rsta
% raden sen n�sta rad tills den har kommit till v�rdet i.
f=figure(402);
n=1;
for i = 1:length(spalt_N)
    for j = 1:length(space)
        subplot(length(spalt_N), length(space), n)
        text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
        plot(abs(fft(spalt(:,i,j))))
        n=n+1;
    end
end

p=1;
for i = 1:length(spalt_N)
    for j = 1:length(space)
        ax = findobj(f,'Type','Axes');
        title(ax(p),['Antal spalter: ' num2str(spalt_N(i)) ' space = ', num2str(space(j))]) % innan s� va [] ist�llet {}, d� br�t den varje g�ng num2str eller ny '' kom.
        ylabel(ax(p),{'Amplitud'})
        xlabel(ax(p),['X ', num2str(0.2*(spalt_N(i)+space(j)-1)/4), ' [cm]' ])
        xlim(ax(p),[-1 (size(spalt(:,i,j),1)-1)/2+1])
        p=p+1;
    end
end
%test = ['L�ngd fr�n spalt till sk�rm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
test = ['Fourie-transform f�r N, spalt\_N och spaces'];
suptitle(test)
%Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures