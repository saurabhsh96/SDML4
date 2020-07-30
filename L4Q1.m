%Lecture 4, Q1
clear;
close all;

%% Defining inputs

% dim
h = 15e-3;
%hs = 2.1e-3;
%er = 12;

% EM
freq = 9e9:0.1e9:11e9;
c = 3e8;
lam = c./freq;
k0 = 2.*pi./lam;

%Additions for Q2
er = 3.46;
hs = lam./(4.*sqrt(er));

% Antenna Dim
L = lam./2; %Halfwavelength
W = lam./20;

%Wave impedances
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;
zeta0 = (sqrt(mu_0/(eps_0*1)));
zetaS = zeta0./(sqrt(er));

%Meshgrid
drad = pi/180;
th = linspace(eps, pi/2-drad, 90);
ph = linspace(eps, 2*pi-drad, 360);
[thi, phi] = meshgrid(th, ph);
dth = thi(1, 2) - thi(1, 1);
dph = phi(2, 1) - phi(1, 1);

%% Q1-1 Propagation constants for TM0, TE1, TM1

%Initial kRho
freqL = length(freq); %largetst freq
kroIm_vec = (-0.2:0.001:0).*k0(end);
kroRe_vec = [(0:0.001:0.5).*k0(end); (0.5:0.001:1).*k0(end)];
%kroRe_vec1 = (0.5:0.001:1).*k0(end);

for ind=1:size(kroRe_vec, 1)
%     [kroRe(ind, :, :), KroIm(ind, :, :)] = meshgrid(kroRe_vec(ind, :), kroIm_vec);
%     kRho(ind, :, :) = kroRe(ind, :, :) + 1j*KroIm(ind, :, :);
%     if(ind == 1)
%        [DlargestTM(ind, :, :), kRReqTM(ind)] = kRhoReq(k0(end), er, h, hs,...
%            squeeze(kRho(ind, :, :)), zeta0, "TM", 1);
%        [DlargestTE, kRReqTE] = kRhoReq(k0(end), er, h, hs,...
%            squeeze(kRho(ind, :, :)), zeta0, "TE", 1);
%     else
%        [DlargestTM(ind, :, :), kRReqTM(ind)] = kRhoReq(k0(end), er, h, hs,...
%            squeeze(kRho(ind, :, :)), zeta0, "TM", 0);
%     end
%Changes for Q2 below
[kroRe(ind, :, :), KroIm(ind, :, :)] = meshgrid(kroRe_vec(ind, :), kroIm_vec);
    kRho(ind, :, :) = kroRe(ind, :, :) + 1j*KroIm(ind, :, :);
    if(ind == 1)
       [DlargestTM(ind, :, :), kRReqTM(ind)] = kRhoReq(k0(end), er, h, hs(end),...
           squeeze(kRho(ind, :, :)), zeta0, "TM", 1);
       [DlargestTE, kRReqTE] = kRhoReq(k0(end), er, h, hs(end),...
           squeeze(kRho(ind, :, :)), zeta0, "TE", 1);
    else
       [DlargestTM(ind, :, :), kRReqTM(ind)] = kRhoReq(k0(end), er, h, hs(end),...
           squeeze(kRho(ind, :, :)), zeta0, "TM", 0);
    end
end

%NR Method
kgPrevTM0 = kRReqTM(2);
kgPrevTM1 = kRReqTM(1);
kgPrevTE = kRReqTE;

%Defining required kgRho
kgRTM0 = zeros(size(freq));
kgRTM1 = zeros(size(freq));
kgRTE = zeros(size(freq));

%Analytical
% kgRTM0A = zeros(size(freq));
% kgRTM1A = zeros(size(freq));
% kgRTEA = zeros(size(freq));

%Intital defination
kgRTM0(freqL) = kRReqTM(2);
kgRTM1(freqL) = kRReqTM(1);
kgRTE(freqL) = kRReqTE;

for ind = size(freq, 2):-1:2
    %Normalized kg
    kgnTM0 = kgPrevTM0./k0(ind);
    kgnTM1 = kgPrevTM1./k0(ind);    
    kgnTE = kgPrevTE./k0(ind);
    
    %Kg for next step
    kgTM0 = kgnTM0.*k0(ind-1);
    kgTM1 = kgnTM1.*k0(ind-1);
    kgTE = kgnTE.*k0(ind-1);
    
    %Required kg
    
    %Changes for Q2
    kgRTM0(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs(ind-1), kgTM0, zeta0, "TM");
    kgRTM1(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs(ind-1), kgTM1, zeta0, "TM");
    kgRTE(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs(ind-1), kgTE, zeta0, "TE");
    
%     kgRTM0(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs, kgTM0, zeta0, "TM");
%     kgRTM1(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs, kgTM1, zeta0, "TM");
%     kgRTE(ind - 1) = findpropSuperStrate(k0(ind-1), er, h, hs, kgTE, zeta0, "TE");
%     
    %Next guess
    kgPrevTM0 = kgRTM0(ind - 1);
    kgPrevTM1 = kgRTM1(ind - 1);
    kgPrevTE = kgRTE(ind - 1);    
end

%Analytical
%kroTE1, kroTM1, kroTM0
[kgRTEA, kgRTM1A, kgRTM0A] =  krhoAnalyticalSuperStrate(k0,lam,h,er);

%Plotting 
%Real part
figure(); %Re(KRho) vs. frequency
plot(freq./10^9, real(kgRTM1A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM1-Analytical'); hold on;
plot(freq./10^9, real(kgRTM1./k0), 'LineWidth', 1.5, 'DisplayName', 'TM1-Numerical');
plot(freq./10^9, real(kgRTE./k0), 'LineWidth', 1.5, 'DisplayName', 'TE-Numerical'); hold on;
plot(freq./10^9, real(kgRTEA./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TE-Analytical');
title(['Guess pole point (Real part) vs. Frequency at \epsilon_r = ', num2str(er)]);
xlabel('Frequency (in GHz)');
ylabel('Re(k_\rho)/k_0');
legend show;
grid on;
hold off;

figure();
plot(freq./10^9, real(kgRTM0A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM0-Analytical'); hold on;
plot(freq./10^9, real(kgRTM0./k0), 'LineWidth', 1.5, 'DisplayName', 'TM0-Numerical');
title(['Guess pole point (Real part) vs. Frequency at \epsilon_r = ', num2str(er)]);
xlabel('Frequency (in GHz)');
ylabel('Re(k_\rho)/k_0');
legend show;
grid on;
hold off;

%Imaginary part
figure(); %Re(KRho) vs. frequency
plot(freq./10^9, imag(kgRTM1A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM1-Analytical'); hold on;
plot(freq./10^9, imag(kgRTM1./k0), 'LineWidth', 1.5, 'DisplayName', 'TM1-Numerical');
plot(freq./10^9, imag(kgRTE./k0), 'LineWidth', 1.5, 'DisplayName', 'TE-Numerical'); hold on;
plot(freq./10^9, imag(kgRTEA./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TE-Analytical');
title(['Guess pole point (Imaginary part) vs. Frequency at \epsilon_r = ', num2str(er)]);
xlabel('Frequency (in GHz)');
ylabel('Im(k_\rho)/k_0');
legend show;
grid on;
hold off;

figure();
plot(freq./10^9, imag(kgRTM0A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM0-Analytical'); hold on;
plot(freq./10^9, imag(kgRTM0./k0), 'LineWidth', 1.5, 'DisplayName', 'TM0-Numerical');
title(['Guess pole point (Imaginary part) vs. Frequency at \epsilon_r = ', num2str(er)]);
xlabel('Frequency (in GHz)');
ylabel('Im(k_\rho)/k_0');
legend show;
grid on;
hold off;

%% 1-2 FF

%Current
M = [1, 0, 0];

%Observation point
r = 1;
z = r.*cos(thi);
z_dash = 0;

%%
%Allocating memory
Eth = zeros([size(k0, 2) size(thi)]);
Eph =  zeros([size(k0, 2) size(thi)]);
Emag =  zeros([size(k0, 2) size(thi)]);
Emax = zeros(size(k0));

%Plotting vec
th_vec = [-thi(1, size(thi, 2):-1:1), thi(1, :)].*180/pi;

for ind = 1:size(k0, 2)
    %Propagation constants
    ks = k0(ind).*sqrt(er);
    kxs = k0(ind).*sin(thi).*cos(phi);
    kys = k0(ind).*sin(thi).*sin(phi);
    kzs = k0(ind).*cos(thi);
    kRho = sqrt(kxs.^2 + kys.^2); 

    %Tx Line Equivalence
    [vTM, vTE, iTM, iTE] = trxline_SuperStrate(k0(ind), er, h, hs,...
        zeta0, zetaS, kRho, z);
    
    %JFT of the current
    %keq = (k0(ind) + ks)./2;
    MFT = CurrentFT(k0(ind), kxs, kys, L(ind), W(ind), M);

    %Calling SGF
    [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(k0(ind), 1, kxs, kys, vTM, ...
        vTE, iTM, iTE, zeta0, zetaS);

    %Calling Field function
    [Eth(ind,:,:), Eph(ind,:,:), Emag(ind,:,:), Emax(ind)] = Field(k0(ind), ...
        kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

    %Free space Field
%     [EFx, EFy, EFz] = FF(freq(ind), 1, L(ind), W(ind), r, thi, phi);
    
    %ErFS = (EFx.*sin(thi).*cos(phi)) + (EFy.*sin(thi).*sin(phi)) + (EFz.*cos(thi));
%     EthFS(ind, :, :) = (EFx.*cos(thi).*cos(phi)) + ...
%         (EFy.*cos(thi).*sin(phi)) - (EFz.*sin(thi));
%     EphFS(ind, :, :) = (-EFx.*sin(phi)) + (EFy.*cos(phi));
%     EmagFS(ind, :, :) = sqrt((abs(EFx)).^2 + (abs(EFy)).^2 + (abs(EFz)).^2); %Magnitude of E
%     EmaxFS(ind) = max(max(EmagFS(ind,:,:)));
%     EmagNormFS = squeeze((EmagFS(ind,:,:)./EmaxFS(ind)));
%     EVecFS = [EmagNormFS(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1), EmagNormFS(1, :)];
%     EVec1FS = [EmagNormFS(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1), EmagNormFS(91, :)];
        
    figure();
    EmagNorm = squeeze((Emag(ind,:,:)./Emax(ind)));
    EVec = [EmagNorm(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1) EmagNorm(1, :)];
    plot(th_vec, mag2db(EVec), 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
    EVec1 = [EmagNorm(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1) EmagNorm(91, :)];
    plot(th_vec, mag2db(EVec1), 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
    title(['Normalized Far-Field vs. \theta, Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    xlabel('\theta (deg)');
    ylabel('Normalized Far-Field E(\theta, \phi) (dB)');
    ylim([-35, 0]);
    legend show;
    hold off;
    grid on;
    
    figure(); %Comparison 0 deg
%    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], [EmagNormFS(1,size(EmagNormFS,2):-1:1), ...
%         EmagNormFS(1,:)], 'LineWidth', 1.5, 'DisplayName', 'Free Space');
%     hold on;
    %polarplot(EVec, 'LineWidth', 1.5, 'DisplayName', 'Spectral Field'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec...
       , 'LineWidth', 1.5, 'DisplayName', 'Spectral');
    title(['Comparison FF FS and Spectral \phi= 0 deg Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
    
    figure(); %Comparison 90 deg
%     polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], [EmagNormFS(91,size(EmagNormFS,2):-1:1), ...
%         EmagNormFS(91,:)], 'LineWidth', 1.5, 'DisplayName', 'Free Space');
%     hold on;
    %polarplot(EVec, 'LineWidth', 1.5, 'DisplayName', 'Spectral Field'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec1...
       , 'LineWidth', 1.5, 'DisplayName', 'Spectral');
    title(['Comparison =FF FS and Spectral \phi= 90 deg Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
    

    figure(); %FF
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec, 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
    polarplot([-thi(1,size(thi,2):-1:1),thi(1,:)], EVec1, 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
    title(['Normalized Far-Field, Freq = ', num2str(freq(ind)./10^9), ' GHz']);
    legend show;
    hold off;
end

%% Q1-3 Propagation constant changing w.r.t er ? Numerical in the same way as previous problem?

er = 2:0.2:25;
freq = 10e9;
lam = c./freq;
k0 = 2.*pi./lam;
hs = lam./(4.*sqrt(er));

%Analytical
[kgRTEA, kgRTM1A, ~] = krhoAnalyticalSuperStrate(k0,lam,h,er);

%Numerical
kroIm_vec = (-0.2:0.001:0).*k0;
kroRe_vec = (0:0.001:0.5).*k0; 

[kroRe, KroIm] = meshgrid(kroRe_vec, kroIm_vec);
kRho = kroRe + 1j*KroIm;

[DlargestTM, kRReqTM] = kRhoReq(k0, er(end), h, hs(end), kRho, zeta0, "TM", 1);
[DlargestTE, kRReqTE] = kRhoReq(k0, er(end), h, hs(end), kRho, zeta0, "TE", 1);

kgRTE = zeros(size(er));
kgRTM = zeros(size(er));

kgPrevTM = kRReqTM;
kgPrevTE = kRReqTE;

kgRTE(end) = kgPrevTE; 
kgRTM(end) = kgPrevTM;

for ind = size(er, 2):-1:2
    kgnTM = kgPrevTM./k0;    
    kgnTE = kgPrevTE./k0;
    
    %Kg for next step
    kgTM = kgnTM.*k0;
    kgTE = kgnTE.*k0;
    
    %Required kg
    kgRTM(ind - 1) = findpropSuperStrate(k0, er(ind - 1), h, hs(ind - 1), kgTM, zeta0, "TM");
    kgRTE(ind - 1) = findpropSuperStrate(k0, er(ind - 1), h, hs(ind - 1), kgTE, zeta0, "TE");
    
    %Next guess
    kgPrevTM = kgRTM(ind - 1);
    kgPrevTE = kgRTE(ind - 1);
end

figure(); %Re(KRho) vs. frequency
plot(er, real(kgRTM1A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM1-Analytical'); hold on;
plot(er, real(kgRTEA./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TE-Analytical');
plot(er, real(kgRTE./k0), 'LineWidth', 1.5, 'DisplayName', 'TE-Numerical');
plot(er, real(kgRTM./k0), 'LineWidth', 1.5, 'DisplayName', 'TM-Numerical');
title(['Guess pole point (Real part) vs. \epsilon_r at Freq = ', num2str(freq./10^9), ' GHz']);
xlabel('\epsilon_r');
ylabel('Re(k_\rho)/k_0');
xlim([2, 25]);
legend show;
grid on;
hold off;

%Imaginary part
figure(); %Re(KRho) vs. frequency
plot(er, imag(kgRTM1A./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TM1-Analytical'); hold on;
plot(er, imag(kgRTEA./k0), '--', 'LineWidth', 1.5, 'DisplayName', 'TE-Analytical');
plot(er, imag(kgRTE./k0), 'LineWidth', 1.5, 'DisplayName', 'TE-Numerical');
plot(er, imag(kgRTM./k0), 'LineWidth', 1.5, 'DisplayName', 'TM-Numerical');
title(['Guess pole point (Imaginary part) vs. \epsilon_r at Freq = ', num2str(freq./10^9), ' GHz']);
xlabel('\epsilon_r');
ylabel('Im(k_\rho)/k_0');
legend show;
xlim([2, 25]);
grid on;
hold off;

%% 1-4 BW of above problem

%er-freq meshgrid
freq = 5e9:0.01e9:15e9;
er = 2:0.05:25;

%[freq, er] = meshgrid(freq_vec, er_vec);
lam = c./freq;
k0 = 2.*pi./lam;

% Antenna Dim
L = lam./2; %Halfwavelength
W = lam./20;

BW = zeros(size(er));
Dm = zeros([size(er, 2) size(freq, 2)]);
% figure();
for ind = 1:length(er)
    for indF = 1:length(freq)
        ks = k0(indF).*sqrt(er);
        kxs = k0(indF).*sin(thi).*cos(phi);
        kys = k0(indF).*sin(thi).*sin(phi);
        kzs = k0(indF).*cos(thi);
        kRho = sqrt(kxs.^2 + kys.^2); 
        
        hs = lam(indF)./(4.*sqrt(er(ind)));
        zetaS = zeta0./(sqrt(er(ind)));
        
        %Tx Line Equivalence
        [vTM, vTE, iTM, iTE] = trxline_SuperStrate(k0(indF), er(ind), h, hs,...
            zeta0, zetaS, kRho, z);

        %JFT of the current
        %keq = (k0(ind) + ks)./2;
        MFT = CurrentFT(k0(indF), kxs, kys, L(indF), W(indF), M);

        %Calling SGF
        [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(k0(indF), 1, kxs, kys, vTM, ...
            vTE, iTM, iTE, zeta0, zetaS);

        %Calling Field function
        [~, ~, Emag, Emax] = Field(k0(indF), ...
            kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);
        
        %Directivity
        [Dir, ~] = DirectivityF(Emag, 1, r, thi, dth, dph);
        
        Dm(ind, indF) = max(max(Dir(:,1)));
    end
    dirDb = pow2db(Dm(ind, :));
    distance = dirDb(2) - dirDb(1);
    Max = max(dirDb);
    dir3Db = Max - 3;
    index = find(ismembertol(dirDb, dir3Db, distance));
    
    %LowerFreq
    Fl = freq(index(1));
    
    %HigherFreq
    Fh = freq(index(end));
    
    %BW
    BW(ind) = 200*(Fh - Fl)./(Fh + Fl);
end

%%
figure(1);
for ind = 1:length(er)
    dirDb = pow2db(Dm(ind, :));
    distance = dirDb(2) - dirDb(1);
    Max = max(dirDb);
    dir3Db = Max - 3;
    index = find(ismembertol(dirDb, dir3Db, 4*distance));
    
    %LowerFreq
    Fl = freq(index(1));
    
    %HigherFreq
    Fh = freq(index(end));
    
    %BW
    BW(ind) = 200*(Fh - Fl)./(Fh + Fl);
    if(mod(er(ind), 4)==0 || er(ind) == 25)
        plot(freq./10^9, dirDb, 'LineWidth', 1.5,'DisplayName', ['\epsilon_r = ', num2str(er(ind))]); hold on;
    end
end
hold off;
grid on;
xlim([9, 11]);
title('Directivity vs. Frequency at different \epsilon_r');
legend show;
xlabel('Frequency (in GHz)');
ylabel('Directivity (in dB)');


figure(2)
This = load('FinalFinal.mat');
B1 = smoothdata(This.BW);
B = smoothdata(BW);
plot(This.er, B1, 'LineWidth', 1.5); hold on;
plot(er, B, 'LineWidth', 1.5);
xlim([4, 25]);
title('Relative BW vs. \epsilon_r');
ylabel('BW(%)');
xlabel('\epsilon_r');
grid on;

%% 