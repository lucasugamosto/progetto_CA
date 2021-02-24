%Luca Sugamosto, matricola 0252792
%Progetto di CA della parte a tempo discreto

clc
clear
close all

s = tf([1 0],1);

%Definizione dell'impianto nel dominio di Laplace
Ps = (0.1*(s+3))/(s*(5*s+1));

%Calcolo i poli della funzione di trasferimento P(s)
[num,den] = tfdata(Ps);
disp("Poli della funzione d'impianto P(s) :");
roots(den{1})

%Inizializzazione del tempo di campionamento Ts, e del tempo di assestamento
%Ta che si vuole ottenere per il sistema a ciclo chiuso finale
Ta = 0.5;                   %Ta = 0.5 s
Ts = Ta/30;                 %Ts = 0.01667 s, (30 = n° campioni per secondo)

%Calcolo della funzione di trasferimento equivalente a tempo discreto P(z)
%dell'impianto da controllare
Pz = c2d(Ps,Ts,'zoh');      %H*P(z)

%Calcolo della funzione di trasferimento P(w) ottenuta da P(z), questo
%grazie al fatto che z e w sono legate dalla trasformazione bilineare
%z = (1+w)/(1-w) e w = (z-1)/(z+1)
Pw = d2c(Pz,'tustin');      %H*P(w)

%Controllo i poli della funzione di trasferimento P(w)
[num,den] = tfdata(Pw);
disp("Poli della P(w): ");
roots(den{1})

%Andando a controllare i poli della funzione P(w), sulla riga di comando
%vengono stampati gli stessi della funzione P(s) e quindi si vede che tale
%funzione ha un polo in 0.
%In realtà andando a studiare il vettore "den" questo ha un polo in
%-4*10^(-13) quindi non è proprio 0. 
%Nonostante ciò, considero quel polo in 0 e quindi non ne devo aggiungere
%un altro nel controllore per asservire segnali costanti di riferimento

%%
%Con questo diagramma di bode vado a vedere fino a quale valore della
%frequenza, la funzione di trasferimento P(s) e P(w) si comportano in modo
%simile e quindi posso dire che fino a quel punto l'approssimazione è
%precisa. Dopo tale valore di frequenza le due funzioni si comportano in
%modo diverso e quindi si può affermare che P(w) non può essere più
%considerata un'approssimazione ottima della P(s)
figure(1)
bode(Ps,Pz,Pw)
legend("Ps","Pz","Pw")
grid on

%%
%Controllo valori dei margini di guadagno e di fase, le frequenze
%caratteristiche e la posizione degli zeri e dei poli della funzione P(w)
figure(2)
margin(Pw)
legend("P(w)");
grid on

figure(3)
rlocus(Pw)
grid on

%%
%Progettazione del controllore C(w) utilizzando i metodi usati anche nel
%caso delle funzioni nel dominio di Laplace
tau = 1/1;                        %zero in -1
alpha = 1/(tau*10);               %Polo in -10  
omega_tau = 1/(tau*sqrt(alpha))   %frequenza a cui si ha lo sfasamento massimo
                                  %e corrisponde al punto medio tra pulsazioni 
                                  %dello zero e del polo
Ra = tf([tau 1],[tau*alpha 1]);

%%
%Vedo cosa succede con questo primo controllo dato solo dalla rete
%anticipatrice
Lw1 = Ra*Pw;
Lw1 = minreal(Lw1);

figure(4)
rlocus(Lw1)
grid on

%%
Cw = 150*Ra;                      %vado a moltiplicare per il guadagno 150
                                  %cosi' da ottenere tutti i poli della
                                  %funzione a ciclo chiuso con parte
                                  %immaginaria nulla
Lw = Cw*Pw; 

%Controllo i poli e gli zeri della funzione d'anello aperto e i possibili poli
%e zeri della funzione a ciclo chiuso
figure(5)
rlocus(Lw,[0 1])
grid on

%Controllo il margine di fase e di guadagno della funzione d'anello aperto
figure(6)
margin(Lw)
legend("L(w)");
grid on

%%
%Calcolo della funzione di sensitività complementare Wyr(w) e calcolo dei
%poli di tale funzione per vedere se sono tutti con parte reale negativa e
%quindi se il sistema a ciclo chiuso è asintoticamente stabile
Wyrw = Lw/(1+Lw);
Wyrw = minreal(Wyrw);

[num,den] = tfdata(Wyrw);
disp("Poli della Wyr(w): ");
roots(den{1})

%Controllo la risposta del sistema a ciclo chiuso Wyr(w) allo step e
%verifico quanto vale il valore del tempo di assestamento Ta
figure(7)
step(Wyrw,5)
legend("Wyr(w)");
grid on
%Nel grafico si vede come è presente una sottoelongazione di valore -0.333
%e si ottiene un tempo di assestamento pari a 0.453 secondi

%Trasformazione del controllore C(w) nel controllore digitale C(z)
Cz = c2d(Cw,Ts,"tustin");

%non serve che vada a controllare Wyr(z) perchè grazie al cambio di variabili
%secondo Tustin la stabilità è mantenuta e quindi se Wyr(w) è
%asintoticamente stabile, lo sarà anche Wyr(z)

%%
%simulazione del controllore e quindi di tutto l'impianto a tempo discreto
%Realizzazione nello spazio di stato del processo a tempo discreto
[A,B,C,D] = ssdata(Pz);
[Ac,Bc,Cc,Dc] = ssdata(Cz);
%A è una matrice 2x2 mentre Ac è una matrice 1x1, da questo definisco il
%numero di righe dei vettori di stato 'x' e 'xc'

N = 250;            %numero di campioni

%inizializzo stato, uscita, ingresso e errore dell'impianto P(z)
x = zeros(2,N);
y = zeros(1,N);
u = zeros(1,N);
e = zeros(1,N);

%inizializzo stato del controllore C(z)
xc = zeros(1,N);

%inizializzo valori delle variabili all'istante iniziale e il riferimento
%che si intende asservire
r = 1;              %sta ad indicare il riferimento a gradino unitario
u(1) = 0;
y(1) = C*x(:,1) + D*u(1);
e(1) = r - y(1);
%x(:,1) = [0;0];

%%
%simulo per istanti di tempo che vanno da 2 a N
for k = 2:N
    %simulazione del processo
    x(:,k) = A*x(:,k-1) + B*u(k-1);
    y(k) = C*x(:,k) + D*u(k-1);
    
    %sistema di acquisizione (e(k) ingresso del controllore)
    e(k) = r - y(k);
    
    %aggiornamento della legge di controllo
    xc(:,k) = Ac*xc(:,k-1) + Bc*e(k-1);
    
    %invio al DAC (digital-analogic converter)
    u(k) = Cc*xc(:,k) + Dc*e(k);
end

%Vado a graficare la variazione nel tempo, dell'uscita del sistema digitale
%per vedere se questa si comporta come la funzione del sistema a ciclo 
%chiuso Wyr(w).
%Se così fosse allora posso affermare di aver progettato un controllore
%digitale che renda il sistema a ciclo chiuso stabile e soddisfa le
%specifiche richieste inizialmente
time = Ts*[0:1:N-1];
figure(8)
plot(time,y,'b','LineWidth',2)
legend("y(k)");
grid on