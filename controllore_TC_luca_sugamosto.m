clc
clear
close all

%Specifiche a regime da soddisfare: ERRORE NULLO PER RIFERIMENTI A RAMPA;
%Specifiche a transitorio da soddisfare: SOVRAELONGAZIONE MINORE DEL 20%
                                        %TEMPO DI ASSESTAMENTO MINORE DI 2s

%Inizializzazione delle variabili
tau = 0.3;           %costante di tempo [s]
b = 50;              %coefficiente di attrito [N*s/m]
m = 1000;            %massa dell'automobile [kg]

s = tf([1 0],1);
%Funzione di trasferimento dell'impianto P(s)
P = 1/(m*tau*s^2 + (b*tau+m)*s + b);

%Controllo i poli della funzione dell'impianto
[num,den] = tfdata(P);
disp("Poli della funzione d'impianto P(s): ");
roots(den{1})

%%
%Progettazione del controllore per avere errore nullo a regime per segnale
%di riferimento a rampa.
%Siccome la funzione P(s) non ha nessun polo nell'origine ne devo mettere
%due al controllore, mentre il guadagno, che è libero, lo pongo unitario.
C1 = 1/(s^2);

%Definizione della funzione di trasferimento d'anello aperto considerando 
%solo il controllore C1(s)
L1 = C1*P;

%Controllo i poli della funzione d'anello aperto
[num,den] = tfdata(L1);
disp("Poli della funzione d'anello aperto L1(s): ");
roots(den{1})

%Controllo l'andamento dei poli e degli zeri all'aumentare del guadagno
%attraverso la funzione rlocus
figure(1)
rlocus(L1)
legend("L1(s)")
grid on
%Da tale grafico si nota che che i due poli posti in zero tendono da subito
%verso il semipiano reale positivo e quindi devo continuare a progettare il
%controllore

%Controllo i valori del margine di fase e di guadagno, la omega di taglio e
%quella critica
figure(2)
margin(L1)
legend("L1(s)")
grid on

%%
%Progettazione del controllore per avere la stabilizzabilità nel sistema a
%ciclo chiuso, definito dalla funzione Wyr(s)
tauRa = 25;                %-> zero della rete si trova in -0.04
alpha1 = 1/(tauRa*30);     %-> polo della rete si trova in -30
alpha2 = 1/(tauRa*40);     %-> polo della rete si trova in -40

%Calcolo i valori di frequenza che rappresentano il punto medio tra le
%pulsazioni associate allo zero e al polo della rete anticipatrice.
%Questi valori coincidono alla frequenza in cui viene applicato lo
%sfasamento massimo della rete anticipatrice.
omega_t_Ra1 = 1/(tauRa*sqrt(alpha1));
omega_t_Ra2 = 1/(tauRa*sqrt(alpha2));

%Progetto le reti anticipatrici che utilizzerò per stabilizzare il sistema
%a ciclo chiuso Wyr(s)
Ra1 = (1 + tauRa*s)/(1 + tauRa*alpha1*s);
Ra2 = (1 + tauRa*s)/(1 + tauRa*alpha2*s);

%Vado ad inserire due reti anticipatrici con stesso zero perchè cosi'
%attirano a se entrambi i poli nulli inseriti con C1(s).
%Mentre vado ad inserire due reti anticipatrici con poli diversi perchè in
%questo modo la parte immaginaria non è subito diversa da zero, ma aumenta
%dopo aver raggiunto un certo guadagno

%Il coefficiente 1.23 lo aggiungo dopo aver effettuato per la prima volta
%la funzione rlocus e aver studiato il guadagno ottimo per avere tutti
%i poli sull'asse reale negativo con parte immaginaria più piccola
C2 = 1.23*Ra1*Ra2*C1;
C2 = minreal(C2);

L2 = C2*P;
L2 = minreal(L2);

[num,den] = tfdata(L2);
disp("poli della funzione d'anello aperto L2(s): ");
roots(den{1})

%Controllo l'andamento dei poli e degli zeri all'aumentare del guadagno
figure(3)
rlocus(L2)
legend("L2(s)")
grid on

%Controllo i valori del margine di fase e di guadagno, la omega di taglio e
%quella critica
figure(4)
margin(L2)
legend("L2(s)")
grid on

%%
%Controllo se la funzione del sistema a ciclo chiuso è as. stabile
%verificando la posizione dei poli della stessa
Wyr = L2/(1 + L2);      
Wyr = minreal(Wyr);

[num,den] = tfdata(Wyr);
disp("poli della funzione del sistema a c.c. Wyr(s): ");
roots(den{1})
%Si ha che tutti i poli tranne due (quelli più vicini all'origine), hanno
%parte immaginaria nulla, mentre gli altri due hanno parte immaginaria non
%nulla ma comunque molto piccola

figure(5)
step(Wyr)
legend("Wyr(s)");
grid on

%%
%Progettazione del filtro di feedforward per soddisfare le specifiche nel
%transitorio e nello specifico il tempo di assestamento, aggiungendo poli
%veloci
p = 1/52;
Ff1 = (1/P)*(1/(((p*s+1)^2)))      %poli doppi veloci in -52
Ff1 = minreal(Ff1);

[num,den] = tfdata(Ff1);
disp("poli della funzione del filtro di feedforward Ff1(s): ");
roots(den{1})

WyrFf1 = ((Ff1*P) + (C2*P))/(1 + (C2*P));
WyrFf1 = minreal(WyrFf1);

[num,den] = tfdata(WyrFf1);
disp("poli della funzione di sensitività complementare WyrFf1(s): ");
roots(den{1})

%Progettazione di un'altro feedforward per verificare la differenza con il
%primo, aumentando in questo il valore dei poli veloci
p = 1/65;
Ff2 = (1/P)*(1/(((p*s+1)^2)))      %poli veloci in -65
Ff2 = minreal(Ff2);

WyrFf2 = ((Ff2*P) + (C2*P))/(1 + (C2*P));
WyrFf2 = minreal(WyrFf2);

[num,den] = tfdata(WyrFf2);
disp("poli della funzione di sensitività complementare WyrFf2(s): ");
roots(den{1})

%%
%Controllo i valori del tempo di assestamento e della sovraelongazione con
%la funzione step applicata alla funzione di sensitività complementare sia
%nel caso semza filtro sia nel caso con filtro
figure(6)
step(Wyr)
hold on
step(WyrFf1)
hold on
step(WyrFf2)
legend("Wyr","WyrFf1","WyrFf2");
grid on

%Per proseguire con i conti ho scelto la funzione di sensitività
%complementare WyrFf1 poichè anche se ha un tempo di assestamento maggiore
%rispetto a WyrFf2, comunque riesce ad asservire meglio il segnale di
%riferimento.

%%
%Visualizzo l'andamento delle varie funzioni che caratterizzano un sistema
%a ciclo chiuso, cioè le funzione WyrFf1, Wer (coincide con Wetruer poichè
%H(s) = 1)
Wer = (1 - (P*Ff1))/(1 + (P*C2));
Wer = minreal(Wer);

W = WyrFf1 + Wer;
W = minreal(W);

figure(7)
bode(WyrFf1,Wer,W)
legend("WyrFf1","Wer","WyrFf1+Wer");
grid on

%%
%Simulazione della risposta del sistema a ciclo chiuso a un segnale di
%riferimento a rampa R/s^2
figure(8)
time = [0:0.01:75];
R = 1;
u = R*time;
lsim(Wyr,u,time)
hold on
lsim(WyrFf1,u,time)

legend("Wyr","WyrFf1")
grid on

%%
%Calcolo del massimo ritardo ammissibile
Mf = (72.5*pi)/180;          %margine di fase di L2(s) [rad]
omega_tau = 0.75;            %frequenza di taglio di L2(s) [rad/s]
disp("massimo ritardo ammissibile: ");
Rmax = Mf/omega_tau          %massimo ritardo ammissibile [s]

%%
%Verifica della robustezza usando il criterio di Kharitonov con variazione
%del fattore di smorzamento e della costante di tempo  che hanno
%una tolleranza del 10% -> b = [45;55], tau = [0.27;0.33]
b_v = [45;55];
tau_v = [0.27;0.33];

%d = denominatore, n = numeratore della funzione associata
Pn = 1;
Pd1 = (m*tau_v(1))*s^2 + (b_v(1)*tau_v(1)+m)*s + b_v(1);
Pd2 = (m*tau_v(2))*s^2 + (b_v(2)*tau_v(2)+m)*s + b_v(2);

Cn = (9.225*10^5)*s^2 + (7.38*10^4)*s + 1476;
Cd = s^4 + 70*s^3 + 1200*s^2;

den_WyrFf1_1 = (Pd1*Cd)+(Pn*Cn);
den_WyrFf1_2 = (Pd2*Cd)+(Pn*Cn);

%%
%per ogni polinomio che ho definito in maniera diversa, vado a vedere se è
%di Hurwitz, se lo sono tutti e quattro allora anche il denominatore di
%Wyr(s) è di Hurwitz e quindi asintoticamente stabile

%funct_1 = {0+ 1+ 2- 3- ...} primo modo di scelta dei coefficienti
funct_1 = tf([270 2.412*10^4 4.673*10^5 1.218*10^6 976500 73800 1476],1);

disp("applicazione del criterio di Routh-Hurwitz su funct_1: ");
myRouth([270 2.412*10^4 4.673*10^5 1.218*10^6 976500 73800 1476])

%funct_2 = {0- 1- 2+ 3+ ...} secondo modo di scelta dei coefficienti
funct_2 = tf([330 1.991*10^4 3.949*10^5 1.226*10^6 988500 73800 1476],1);

disp("applicazione del criterio di Routh-Hurwitz su funct_2: ");
myRouth([330 1.991*10^4 3.949*10^5 1.226*10^6 988500 73800 1476])

%funct_3 = {0+ 1- 2+ 3- ...} terzo modo di scelta dei coefficienti
funct_3 = tf([330 1.991*10^4 4.673*10^5 1.218*10^6 988500 73800 1476],1);

disp("applicazione del criterio di Routh-Hurwitz su funct_3: ");
myRouth([330 1.991*10^4 4.673*10^5 1.218*10^6 988500 73800 1476])

%funct_4 = {0- 1+ 2- 3+ ...} quarto modo di scelta dei coefficienti
funct_4 = tf([270 2.412*10^4 3.949*10^5 1.226*10^6 976500 73800 1476],1);

disp("applicazione del criterio di Routh-Hurwitz su funct_4: ");
myRouth([270 2.412*10^4 3.949*10^5 1.226*10^6 976500 73800 1476])