clear; 
close all; 
clc;
count = zeros(1,9);
%NMSE = zeros(1,9);
for w = 1:9
for iter = 1:1000 
%% System Params
L       = 2;                    % number of paths (including LOS)
Rs      = 100e6;                % total BW in Hz
N       = 10;                   % number of subcarriers 
Nt      = 16;                   % number of TX antennas
Nr      = Nt;                   % number of RX antennas
Nb      = Nt*2;                 % number of beams in dictionary
Ns      = 5;                    % number of beams sent (randomly beams)
c       = 3e8;                  % speed of light m/s
Ts      = 1/Rs;                 % sampling period in us
posTx   = [0 0]';               % TX is assumed to be in [0, 0]
posRx   = [4 1]';               % RX (user) position
SP      = [2, 2];               % scatter point position
alpha   = 0.2;                  % user orientation
h       = 10*ones(1,L);         % channel gain
GR      = (sqrt(5) - 1)/2;      % Golden number
%sigma   = 0.1;                  % noise std
SNR_db = -10:5:30;
W = length(SNR_db);
%% Compute Channel Parameters for L paths
TOA = zeros(1, L); AOD = zeros(1, L); AOA = zeros(1, L);
TOA(1) = norm(posRx)/c;                                                     % LOS TOA
AOD(1) = atan2(posRx(2) - posTx(2), posRx(1) - posTx(1));                   % LOS AOD
AOA(1) = atan2(posTx(2) - posRx(2), posTx(1) - posRx(1)) - alpha;           % LOS AOA
for p = 1:L-1
    TOA(p+1) = (norm(SP(p,:)) + norm(posRx - SP(p,:)'))/c;                  % NLOS TOA
    AOD(p+1) = atan2(SP(p,2), SP(p,1));                                     % NLOS AOD
    AOA(p+1) = atan2(SP(p,2) - posRx(2), SP(p,1) - posRx(1)) - alpha;       % NLOS AOD
end
%% Create dictionary
Ut = zeros(Nt,Nb);
Ur = zeros(Nr,Nb);
aa = -Nb/2:Nb/2-1;
aa = 2*aa/Nb;
for m = 1:Nb
    Ut(:,m) = getResponse(Nt,aa(m))*sqrt(Nt);
    Ur(:,m) = getResponse(Nr,aa(m))*sqrt(Nr);
end
%% Generate channel
H = zeros(Nr,Nt,N); A_rx = zeros(Nr,L); A_tx = zeros(Nt,L); Gamma = zeros(L, L, N);
Hb = zeros(Nb, Nb, N);
for n = 1:N
    for p = 1:L
        A_rx(:,p) = getResponse(Nr,sin(AOA(p)))*sqrt(Nr);
        A_tx(:,p) = getResponse(Nt,sin(AOD(p)))*sqrt(Nt);
        Gamma(p,p,n) = h(p)*exp(-1j*2*pi*TOA(p)*(n-1)/(N*Ts));
        H(:,:,n) = H(:,:,n) + A_rx(:,p)*Gamma(p,p,n)*A_tx(:,p)';
    end
    Hb(:,:,n) = Ur'*H(:,:,n)*Ut;
end

%% Visualize the beamspace channel for 1 subcarrier in AOA/AOD space
Hb = zeros(Nb, Nb, N);
for n = 1:N
    Hb(:,:,n) = Ur'*H(:,:,n)*Ut;
end
Hb_mat = reshape(Hb, [Nb*Nb  N]);
Hb_vec = reshape(Hb, [Nb*Nb*N  1]);
% figure;
% subplot(1, 2, 1);
Hb=Ur'*H(:,:,1)*Ut;
% mesh(abs(Hb));
% xlabel('AOD'); ylabel('AOA');

%% Generate the observation and beamformers
y = zeros(Nr,Ns,N); 
signal = zeros(Nr,Ns,N); 
noise = zeros(Nr,Ns,N); 
F = zeros(Nt,Ns,N);
for g = 1:Ns    
    for n = 1:N        
        F(:,g,n) = exp(1j*rand(Nt,1)*2*pi);                                % random beamformers (note: we don't add data symbols, they are part of F)
        %F(:,:,n) = F(:,:,n)/norm(F(:,:,n), 'fro'); % normalize power of F to 1 (fro norm)
        signal(:,g,n) = H(:,:,n)*F(:,g,n); 
    end
end
signal_vec = signal(:);

for n = 1:N
    for n = 1:W
    PoS(n) = 1/(Nr*Ns) * sum(abs(signal_vec).^2);    
    PoN(n) = PoS(n)/(10^(SNR_db(n)/10));
    sigma(n) = sqrt(PoN(n));
    end
end

for g = 1:Ns
    for n = 1:N
        for n = 1:W
        noise(:,g,n) = sigma(n)/sqrt(2)*(randn(Nr,1) + 1i*randn(Nr,1));        % noise
        y(:,g,n) = signal(:,g,n) + noise(:,g,n);  
        end
    end
end


%% Vectorize and generation of the basis
yb = zeros(Nr*Ns,N);
Omega = zeros(Nr*Ns,Nb*Nb,N);
for n = 1:N
    yb(:,n) = reshape(y(:,:,n), Nr*Ns,1);
    nb(:,n) = reshape(noise(:,:,n), Nr*Ns,1);
    Omega(:,:,n) = kron((Ut'*F(:,:,n)).',Ur);
end
y_vec = reshape(yb, [Nr*Ns*N, 1]);
noise_vec = reshape(nb, [Nr*Ns*N, 1]);
%test = norm(yb(:,3) - (1/Nb)^2*Omega(:,:,3)*Hb_mat(:,3))

%% DCS-SOMP
[indices,h_hat]=DCSSOMP(yb(:,w),Omega,L);                  % the last input is the number of paths it recovers 
support_set = sort(abs(indices));
%for i=1:L
   % support_set(i) = find(abs(indices) == h(Nb*Nb - L + i));
%end
%support_set = sort(support_set)
if (support_set(1) == 656 && support_set(2) == 891)
    count(w) = count(w) + 1
end
end
hb_est = zeros(Nb*Nb, 1);
%hb_est(Omega_k) = x_k;
for  i = 1:length(indices)
    hb_est(indices(i)) = Hb(indices(i))
end
Hb_est = reshape(hb_est, Nb, Nb);
%H_b = (Hb./(Nt*Nr));
%NMSE(w) = norm(Hb_est - Hb , 'fro')/norm(Hb, 'fro')
end
% Plot
figure(1)
subplot(121)
mesh(asin(aa),asin(aa),abs(Hb(:,:,1)));
title('Channel in Angular Domain');
xlabel('AOD (rad)'); ylabel('AOA (rad)');
subplot(122)
mesh(asin(aa),asin(aa), abs(Hb_est));
title('Estimated channel in Angular Domain by SOMP')
xlabel('AOD (rad)'); ylabel('AOA (rad)'); zlabel('|Gain|');
toc

