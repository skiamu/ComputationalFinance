function [ price, check_marting ] = FFT(Strike, S0, T, r, model, p)
%calcola il prezzo di una European Call option con il metodo di Carr &
%Madan e il modello B&S
%   INPUT:
%         Strike = vettore di strike su cui si vuole calcolare il prezzo
%         S0 = valore iniziale sottostante
%         T = time to maturity
%         r = risk-free interest rate
%         model = stringa che indica il modello 
%         p = vettore coi parametri del modello
%         
%   OUTPUT:
%         price = vettore di prezzi della opzione Call
%         check_marting = controllo sulla condizione di martingalità, deve
%                         essere 1
%
%   uses: Characteristic_Exp

tic
%% setting parameters

Npow = 15;
% numero di punti per la quadratura
N = 2^Npow;
% troncatura dominio di integrazione
A = 1200;
% passo griglia integrazione
eta = A/N;
% discretizzazione dominio integrazione, % non zero altrimenti NaN
% faccio (N-1)/N per avere esattamente N punti e non N+1
v = 0:eta:(N-1)/N*A;   v(1)=1e-22;

lambda =  (2 * pi) / (N * eta);
% log-strike grid (N punti)
k = -lambda * N / 2 + lambda * (0:N-1);

% la funzione caratteristica di Heston non si scrive Psi * T quindi si
% tratta a parte
if strcmp(model,'Heston')
   
   Psi_X = Characteristic_Exp( model, p );
   
   Phi_X = @(u) exp(Psi_X(T,u));
    
else
   
   % esponente caratteristico
   Psi_X_star = CharacteristicExp( model, p );
   
   % funzione caratteristica in cui si è preso come drift b = - Psi_X_star(-i)
   % per soddisfare la condizione di martingalità
   Phi_X = @(v) exp(T * (Psi_X_star(v) - 1i * v * Psi_X_star(-1i)));
   
end

% Phi_X(-1i) = 1 per condizione martingala
check_marting = Phi_X(-1i);

%% FFT application

% trasformata di X_T calcolata nei nodi di integrazione
Z_k = (exp(1i .* v * r * T) ./ (1i * v .* (1 + 1i * v))) .* (Phi_X(v - 1i) - 1);
% pesi per quadratura (trapezi)
w = ones(1,N);   w(1) = 0.5;    w(end) = 0.5;
% integranda FFT calcolata nel v_j
x_k = w .* eta .* Z_k .* exp(1i * pi * (0:N-1));

z_T = (1 / pi) * real(fft(x_k));
% prezzo Call option per unità di sottostante
c_k = z_T + max(1 - exp(k - r * T), 0);
% prezzo Call option
C = S0 * c_k;
% stike per cui si ha il prezzo
K = S0 * exp(k);


%% Output

% considero solo strike in range sensato
idx = find(K>0.1*S0 & K<3*S0);
C = C(idx);
K = K(idx);
price = interp1(K, C, Strike, 'spline');
elapsed = toc;



end

