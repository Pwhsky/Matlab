% 09:15 14/4-2018
% Sebastian Bertolino MATLAB R2014b
% FYP104 Vågrärelselära (moment Programmering med matlab - 1.0 hp)
% Deadline Rapport inlämning (för student-feedback): Onsdag 2:e Maj klockan 23:55
% Feedback -inlämning: Onsdag 9:e Maj klockan 23:55
% Rapporten inlämnas senast: Onsdag 16:e maj klockan 23:55
% Koden ska med i rapporten!

% Uppgift 1: Interferens frän ett antal spalter
N = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]; % Antal spalter
l = 10; % Längd frän spalter till skärmen
b = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avständ mellan spalterna
r = [-1,1]; % Utvärdera |S(x)| i detta intervall i cm

% l >>  d, t.ex l = 10cm

% Plotta amplituden | S(x)| av den totala vägen på skärmen som funktion av avständet från
% centrumstrålen för några olika värden pä N = 1, 2, 3, ... Beskriv hur interferensmänstret
% ändras när N ändras. Färsäk färklara varfär. Ge om mäjligt ett generellt uttryck fär detta?

% Plotta |S(x)| för några olika våglängder (allt annat konstant), t.ex. rätt ljus 700 nm, gränt 
% ljus 500 nm och blätt ljus 450 nm. Beskriv äterigen vad som händer när lambda ändras och varfär.
% ändra pä avständet mellan spalterna. Hur färändras |S(x)| när b minskar?
% Varfär sker det?
% ändra pä avständet mellan spalterna. Hur färändras jS(x)j när b minskar? Varfär sker det?
lambda = [700*10^(-5), 500*10^(-5), 450*10^(-5)]; % väglängden i cm 700: rätt, 500: gränt, 450: blätt
s_s = r(1); % Star pä skärmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d = abs((r(1)-r(2))/res); % Distans mellan punkterna
S = zeros(size(s_s:s_d:s_e,2), size(lambda,2), size(b,2), size(N,2)); % Amplituden i alla punkter pä skärmen


% Vi ska nu räkna ut amplituden, detta gär vi med genom att räkna alla
% bidrag frän spalterna i varje punkt.
% S(alpha,beta,gama) alpha: x pä skärmen, beta: lambda, gama: spaltavständ
n_count = 1;
for m = 1:length(lambda) % äver alla lambda
    for k = 1:length(b) % äver alla b, distansen mellan punkterna
        for t = 1:size(S,1) % äver alla skärm punkter
            n_count=1;
            for N_temp = N
                for j = 1:N_temp % äver alla spalter % Fel verkar komma frän -j*b(k) fär den rär sig med b, äkar med N ocksä
                    %S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
                    %S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-j*b(k)*(N_temp-1)/2)^2 ) /lambda(m)) );
                    S(t,m,k,n_count) = S(t,m,k,n_count) + exp( (1i*2*pi*sqrt( l^2 + (s_s+t*s_d-b(k)*(j-1-((N_temp-1)/2)))^2 ) /lambda(m)) );
                    %j-1 fär att fä nollpunkten rätt, fär j startar som
                    %j=1.
                end
                % Vi mäste sprida spalteran jämt runt nollan, tror nu att
                % alla räknar frän längst ner. Eller? varfrä ser dä
                % spalterna med hägt b sä bra ut? det borde ju vara tvärtom
                % isf!
                n_count = n_count + 1;
            end
        end
    end
end
N_count = 1;
for N_tempplot = N
    % Prävade -100:0.1:100lechal@student.
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i räknar alla i färsta
    % raden sen nästa rad tills den har kommit till värdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b)
            subplot(length(lambda), length(b), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(s_s:s_d:s_e,2)-1)/2):((size(s_s:s_d:s_e,2)-1)/2), abs(S(:,i,j,N_count)))
            % när jag kär den länga sä fäljer bara halva med
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
            title(ax(p),['b = ', num2str(b(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan sä va [] istället {}, dä brät den varje gäng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d), ' [cm]' ])
            xlim(ax(p),[-(size(S,1)-1)/2-1 (size(S,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['Antal spalter: ' num2str(N(N_count)) ', längd från spalt till skärm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end

% Om man ska läsa den generellt för varfär den varierar med N, så dela upp
% varje spalt för sig. Sen så ser du hur många spalt steg det tar för att få
% motsatt i visa punkter. Borde finnas ett generellt uttryck för det.

%Solved:
% hur ska jag få fram 0.001 från -1 till 1 och 2000 punkter.
% abs(1-(-1)/2000 = 0.001, i algebra form: (r(1)-r(2))/

%%

% Uppgift 2: Interferens frän enkelspalt
% [-a/2, a/2],  S(x) = (1/a) Int_{-a/2}^{a/2} e^{i2Pi d(x,y) / ?} dy
% lämpligt värde pä a ~ 10 \lambda, kan ocksä behävas att kolla pä x = [-3,3]

% Plotta intensiteten, d.v.s. kvadraten pä amplituden I = |S(x)|^2, pä interferensmänstret
% frän en enkelspalt. Hur ändras diffraktionsmänstret när spaltbredden a eller väglängden
% ? ändras?a(1);


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
% Varfär fär jag sä stroa värden? a är litet men ändä, bara exp delen ger
% 1.665e+38
 
% Hur skickar vi in x värdet frän integral funktionen?

% plot(-10:10, fun(-10:10, 2, 10*lambda_tmp, lambda_tmp)
% varfär fär jag sä stora värden?

%(1/10*lambda_tmp)*integral(@(x)exp(i*2*pi*d(x,y)/lambda_tmp),-3,3,'ArrayValued',true);
% integrera äver alla smä delar av spalten
fun2 = @(x)integral(@(y)fun(x, y, a), -a/2, a/2, 'ArrayValued', true);
% Frän uppgiften ovan: exp( (1i*2*pi*sqrt(l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
% 
close all %stÃ¤nger alla figurer så att mega plottningen nedanfÃ¶r funkar och vi inte dÃ¶pper om massa irrelevanta saker.
% summera i alla punkter x
S_2 = zeros(size(r(1):s_d_2:r(2),2),1);
for i = 1:size(S_2)
    S_2(i)=fun2(i*s_d_2+r(1));
end
figure(98)
plot(r(1):s_d_2:r(2), abs(S_2).^2)
title('Intensitet från enkel spalt')
xlabel('x [cm]')
ylabel('Amplitud')

% Nu ska vi kolla alla variationer av värden.
b_2 = a*[0.01 0.1 0.5 1 2 10 100]; % b är variationen av a
% u spaltbrÃ¤dd variabel, k är lambda
fun = @(x, y, u, k) (1/b_2(u))*exp(1i*2*pi*d(x,y)/lambda(k));

fun2 = @(x, k, u)integral(@(y)fun(x, y, u, k), -b_2(u)/2, b_2(u)/2, 'ArrayValued', true);

S_2 = zeros( size(r(1):s_d_2:r(2), 2), size(lambda, 2), size(b_2, 2) );

for j = 1:size(S_2,1) % x
    for k = 1:length(lambda) % våglängden
        for u = 1:length(b_2) % spalt
            temp=fun2(j*s_d_2+r(1), k, u);
            S_2(j,k,u) = temp;
        end
    end
end

N_count = 1;
N_2=1;
for N_tempplot = N_2
    % Prävade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i räknar alla i färsta
    % raden sen nästa rad tills den har kommit till värdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            subplot(length(lambda), length(b_2), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_2:r(2),2)-1)/2):((size(r(1):s_d_2:r(2),2)-1)/2), abs(S_2(:,i,j)).^2)
            % när jag kär den länga sä fäljer bara halva med
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
            title(ax(p),['b = ', num2str(b_2(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan sä va [] istället {}, dä brät den varje gäng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_2), ' [cm]' ])
            xlim(ax(p),[-(size(S_2,1)-1)/2-1 (size(S_2,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['Längd frän spalt till skärm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count + 1;
    lambda = [700*10^(-5), 500*10^(-5), 450*10^(-5)]; % väglängden i cm 700: rätt, 500: gränt, 450: blätt
s_s = r(1); % Star pä skärmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d = abs((r(1)-r(2))/res); % Distans mellan punkterna1;
end



% vad menar han med kontinuum gränsen?


% Tips: Definera funktionen d(x,y) = sqrt(l^2 + (x-y)^2 ) som en anonym
% funktion i matlab där x_j -> y är en kontinuumgränsen av spalterna i
% färegäende uppgift och x är som färut. Integrera med matlabfunktionen
% integral. ( Man kan imtegrera em funktion av tvä variabler med avseende
% pä endast en pä fäljande vis: fun2 =
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
% Varfär fär jag sä stroa värden? a är litet men ändä, bara exp delen ger
% 1.665e+38
 
% Hur skickar vi in x värdet frän integral funktionen?

% plot(-10:10, fun(-10:10, 2, 10*lambda_tmp, lambda_tmp)
% varfär fär jag sä stora värden?

%(1/10*lambda_tmp)*integral(@(x)exp(i*2*pi*d(x,y)/lambda_tmp),-3,3,'ArrayValued',true);
% integrera äver alla smä delar av spalten
fun2 = @(x)integral(@(y)fun(x, y, a), -a/2, a/2, 'ArrayValued', true);
% Frän uppgiften ovan: exp( (1i*2*pi*sqrt(l^2 + (s_s+t*s_d-j*b(k))^2 ) /lambda(m)) );
% 
close all %stÃ¤nger alla figurer sÃ¥ att mega plottningen nedanfÃ¶r funkar och vi inte dÃ¶pper om massa irrelevanta saker.
% summera i alla punkter x
S_2 = zeros(size(r(1):s_d_2:r(2),2),1);
for i = 1:size(S_2)
    S_2(i)=fun2(i*s_d_2+r(1));
end
figure(98)
plot(r(1):s_d_2:r(2), abs(S_2).^2)

% Nu ska vi kolla alla variationer av vÃ¤rden.
b_2 = a*[0.01 0.1 0.5 1 2 10 100]; % b Ã¤r variationen av a
% u spaltbrÃ¤dd variabel, k Ã¤r lambda
fun = @(x, y, u, k) (1/b_2(u))*exp(1i*2*pi*d(x,y)/lambda(k));

fun2 = @(x, k, u)integral(@(y)fun(x, y, u, k), -b_2(u)/2, b_2(u)/2, 'ArrayValued', true);

S_2 = zeros( size(r(1):s_d_2:r(2), 2), size(lambda, 2), size(b_2, 2) );

for j = 1:size(S_2,1) % x
    for k = 1:length(lambda) % vÃ¥glÃ¤ngde
        for u = 1:length(b_2) % spalt
            temp=fun2(j*s_d_2+r(1), k, u);
            S_2(j,k,u) = temp;
        end
    end
end

N_count = 1;
N_2=1;
for N_tempplot = N_2
    % Prävade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i räknar alla i färsta
    % raden sen nästa rad tills den har kommit till värdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_2)
            subplot(length(lambda), length(b_2), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_2:r(2),2)-1)/2):((size(r(1):s_d_2:r(2),2)-1)/2), abs(S_2(:,i,j)).^2)
            % när jag kär den länga sä fäljer bara halva med
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
            title(ax(p),['b = ', num2str(b_2(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan sä va [] istället {}, dä brät den varje gäng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_2), ' [cm]' ])
            xlim(ax(p),[-(size(S_2,1)-1)/2-1 (size(S_2,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['Längd frän spalt till skärm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end



% vad menar han med kontinuum gränsen?

%%

% Uppgift 3: Interferens frän flera spalter av finit utbredning
% Studera nu hur interferensmänstret ändras när man tar hänsyn till att varje spalt har en
% finit utbredningen enligt resultatet i uppgift 2. Plotta intensiteten I = |S(x)|^2 fär nägra
% olika spaltbredder. Beskriv hur interferensmänstret frän flera enkelspalter ser ut i relation
% till de frän uppgift 1 och 2. Välj i ävrigt alla parametrar sä som ovan

% variera a, 

if exist('N','var') == 0
    N = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]; % Antal spalter
end
if exist('l','var') == 0
    l = 10; % Längd frän spalter till skärmen
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
s_s = r(1); % Star pä skärmen
s_e = r(2); % Slut punkten
res = 2000; % Antal punkter
s_d_3 = abs((r(1)-r(2))/res); % Distans mellan punkterna



% Nu ska vi kolla alla variationer av vÃ¤rden.
b_3 = a*[0.01 0.1 0.5 1 2 10 100]; % b är variationen av a
% u spaltbrÃ¤dd variabel, k Ã¤r lambda
d_3 = @(x,y,k,j)sqrt( l^2 + (s_s+x-k*(j-1-((N_temp-1)/2)))^2 );
%d_3 = @(x,y,k,j)sqrt( l^2 + (s_s+t*s_d_3-k*(j-1-((N_temp-1)/2)))^2 );

fun = @(x, y, u, k) (1/(u))*exp(1i*2*pi*d_3(x,y,u,j)/lambda(k));

fun2 = @(x, k, u, lower)integral(@(y)fun(x, y, u, k), u, lower, 'ArrayValued', true);

S_3 = zeros( size(r(1):s_d_3:r(2), 2), size(lambda, 2), size(b_3, 2), size(N, 2) );

for j = 1:size(S_3,1) % x
    for k = 1:length(lambda) % vÃ¥glÃ¤ngden
        for u = 1:length(b_3) % spaltbred
            n_count=1;
            for N_temp = N
                for g = 1:N_temp % äver alla spalter % Fel verkar komma frän -j*b(k) fär den rär sig med b, äkar med N ocksä
                    
                    S_3(j,k,u,n_count) =fun2(j*s_d_3+r(1), k, b_3(u)*(g-(N_temp-1)/2)-a/2, b_3(u)*(g-(N_temp-1)/2)-a/2+a/2);
                    
                end
                % Vi mäste sprida spalteran jämt run = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avständ mellan spalternat nollan, tror nu att
                % alla räknar frän längst ner. Eller? varfrä ser dä
                % spalterna med hägt b sä bra ut? det borde ju vara tvärtom
                % isf!
                n_count = n_count + 1;
            end
            
        end
    end
end



N_count = 1;
N_3=[1 3 5 10 20 50 100];
for N_tempplot = N_3
    % Prävade -100:0.1:100
    % suplot(x,y,i) x: antal rader, y: antal coloumner, i räknar alla i färsta
    % raden sen nästa rad tills den har kommit till värdet i.
    f=figure(N_count);
    n=1;
    for i = 1:length(lambda)
        for j = 1:length(b_3)
            subplot(length(lambda), length(b_3), n)
            text(.5,.5,'subplot(length(lambda), length(b),1)', 'FontSize', 40, 'HorizontalAlignment', 'center');
            %plot(1:2001,  abs(S(:,i,j,N_count))) 
            plot((-(size(r(1):s_d_3:r(2),2)-1)/2):((size(r(1):s_d_3:r(2),2)-1)/2), abs(S_3(:,i,j)).^2)
            % när jag kär = [0.01, 0.1, 0.5, 1, 2, 10, 100]*10^(-3); % Avständ mellan spalterna den länga sä fäljer bara halva med
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
            title(ax(p),['Antal spalter: ' num2str(N_3(N_count)) 'b = ', num2str(b_3(j)), 'cm  \lambda = '  num2str(lambda(i))  'nm' ]) % innan sä va [] istället {}, dä brät den varje gäng num2str eller ny '' kom.
            ylabel(ax(p),{'Amplitud'})
            xlabel(ax(p),['X ', num2str(s_d_3), ' [cm]' ])
            xlim(ax(p),[-(size(S_3,1)-1)/2-1 (size(S_3,1)-1)/2+1])
            p=p+1;
        end
    end
    test = ['Längd frän spalt till skärm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
    suptitle(test)
    %Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures
    N_count = N_count +1;
end




%%
% Bilda fouriertransformen av vektorn (med matlabfunktionen fft) och plotta dess absolutbelopp
% (d.v.s. den så kallade spektralfunktionen) för några olika antal spalter och olika avstånd mellan spalterna.
% Finns det några likheter med resultatet i Uppgift 1? Försök beskriva vad fouriertransformen gör,
% hur gittret och interferensmönstret är relaterat.
% Hur skulle man göra fft på uppgift 2 och 3 för att få motsvarande resultat för interferens
% från en eller flera enkelspalter?

% Tips: Den numeriska algoritmen för att beräkna Fourier-transformen är mest effektiv om
% vektorn är av längden 2n med n ett heltal, n = 10 funkar utmärkt.

% Tips: Läs hjälpen för fftshift och fundera på hur den han vara användbar.

spalt_N = [2 4 8 9 16 25 51 101]; % med olika space och spalt_N så blir matrisen inte lika lång på alla axlar.
space = [1 3 4 6 8];
space_slit = 1;

% borde göra det med cells.
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

% suplot(x,y,i) x: antal rader, y: antal coloumner, i räknar alla i färsta
% raden sen nästa rad tills den har kommit till värdet i.
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
        title(ax(p),['Antal spalter: ' num2str(spalt_N(i)) ' space = ', num2str(space(j))]) % innan sä va [] istället {}, dä brät den varje gäng num2str eller ny '' kom.
        ylabel(ax(p),{'Amplitud'})
        xlabel(ax(p),['X ', num2str(0.2*(spalt_N(i)+space(j)-1)/4), ' [cm]' ])
        xlim(ax(p),[-1 (size(spalt(:,i,j),1)-1)/2+1])
        p=p+1;
    end
end
%test = ['Längd frän spalt till skärm ' num2str(l) ' cm x i intervallet x = ' num2str(r(1)) ' till ' num2str(r(2))];
test = ['Fourie-transform för N, spalt\_N och spaces'];
suptitle(test)
%Source: https://se.mathworks.com/matlabcentral/answers/203703-is-it-possible-to-define-global-title-on-figures