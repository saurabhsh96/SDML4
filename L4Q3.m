%Lecture 4, Q3
clear;
close all;

%% Input defination

% dim
h = 15e-3;
%hs = 2.1e-3;
er = 12;

% EM
freq = 10e9;
c = 3e8;
lam = c./freq;
k0 = 2.*pi./lam;

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
ph = (0:5:360)*pi/180;
th = linspace(eps,89,181)*pi/180;

dth = th(2)-th(1);
dph = ph(2)-ph(1);

[thi,phi] = meshgrid(th,ph);

M = [1, 0, 0];
r = 1;
z = r.*cos(thi);
z_dash = 0;
%% Q3-1

d = lam./2;

ks = k0.*sqrt(er);
kxs = ks.*sin(thi).*cos(phi);
kys = ks.*sin(thi).*sin(phi);
kzs = ks.*cos(thi);
kRho = sqrt(kxs.^2 + kys.^2); 

%Tx Line Equivalence
[vTM, vTE, iTM, iTE] = trxline_SuperStrate3(k0, er, h,...
    zeta0, zetaS, kRho, z);

%JFT of the current
%keq = (k0(ind) + ks)./2;
MFT = twoAntFT(k0, kxs, kys, L, W, M, d);

%Calling SGF
[Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(ks, er, kxs, kys, vTM, ...
    vTE, iTM, iTE, zeta0, zetaS);

%Calling Field function
[Eth, Eph, Emag, Emax] = Field(ks, ...
    kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);

figure();
th_vec = [-thi(1, size(thi, 2):-1:1), thi(1, :)].*180/pi;
EmagNorm = squeeze((Emag./Emax));
EVec = [EmagNorm(round(size(phi, 1)./2 + 1), size(thi, 2):-1:1) EmagNorm(1, :)];
plot(th_vec, mag2db(EVec), 'LineWidth', 1.5, 'DisplayName', '\phi = 0 deg'); hold on;
EVec1 = [EmagNorm(round(size(phi, 1)./4 + 1), size(thi, 2):-1:1) EmagNorm(91, :)];
plot(th_vec, mag2db(EVec1), 'LineWidth', 1.5, 'DisplayName', '\phi = 90 deg');
title(['Normalized Far-Field vs. \theta, Freq = ', num2str(freq./10^9), ' GHz']);
xlabel('\theta (deg)');
ylabel('Normalized Far-Field E(\theta, \phi) (dB)');
ylim([-50, 0]);
legend show;
hold off;
grid on;

%% Q3-2 BW of antenna

%er-freq meshgrid
freq = 5e9:0.02e9:15e9;
er = 2:0.05:25;

%[freq, er] = meshgrid(freq_vec, er_vec);
lam = c./freq;
k0 = 2.*pi./lam;

% Antenna Dim
L = lam./2; %Halfwavelength
W = lam./20;

Dm = zeros([size(er, 2) size(freq, 2)]);
% figure();
for ind = 1:length(er)
    for indF = 1:length(freq)
        ks = k0(indF).*sqrt(er(ind));
        kxs = ks.*sin(thi).*cos(phi);
        kys = ks.*sin(thi).*sin(phi);
        kzs = ks.*cos(thi);
        kRho = sqrt(kxs.^2 + kys.^2); 
        
        %hs = lam(indF)./(4.*sqrt(er(ind)));
        zetaS = zeta0./(sqrt(er(ind)));
        
        %Tx Line Equivalence
        [vTM, vTE, iTM, iTE] = trxline_SuperStrate3(k0(indF), er(ind), h,...
            zeta0, zetaS, kRho, z);

        %JFT of the current
        d = lam(indF)./2;
        %keq = (k0(ind) + ks)./2;
        MFT = twoAntFT(k0(indF), kxs, kys, L(indF), W(indF), M, d);

        %Calling SGF
        [Gxx, Gyx, Gzx, Gxy, Gyy, Gzy] = SpectralGFemF(ks, er(ind), kxs, kys, vTM, ...
            vTE, iTM, iTE, zeta0, zetaS);

        %Calling Field function
        [~, ~, Emag, Emax] = Field(ks, ...
            kzs, r, thi, phi, Gxx, Gyx, Gzx, Gxy, Gyy, Gzy, MFT, z, z_dash);
        
        %Directivity
        [Dir, ~] = DirectivityF(Emag, er(ind), r, thi, dth, dph);
        
        Dm(ind, indF) = max(max(Dir(:,1)));
    end
%     dirDb = pow2db(Dm(ind, :));
%     distance = dirDb(2) - dirDb(1);
%     Max = max(dirDb);
%     dir3Db = Max - 3;
%     index = find(ismembertol(dirDb, dir3Db, distance));
%     
%     %LowerFreq
%     Fl = freq(index(1));
%     
%     %HigherFreq
%     Fh = freq(index(end));
%     
%     %BW
%     BW(ind) = 200*(Fh - Fl)./(Fh + Fl);
end

%%

for ind = 1:length(er)
    dirDb = pow2db(Dm(ind, :));
    distance = dirDb(2) - dirDb(1);
    Max = max(dirDb);
    dir3Db = Max - 3;
    index = find(ismembertol(dirDb, dir3Db, 5*distance));
    
    %LowerFreq
    Fl = freq(index(1));
    
    %HigherFreq
    Fh = freq(index(end));
    
    %BW
    BW(ind) = 200*(Fh - Fl)./(Fh + Fl);
end

figure();
% prev = load('prev.mat');
B = smoothdata(BW);
% B2 = smoothdata(prev.B);
plot(er, B, 'LineWidth', 1.5, 'DisplayName', 'Substrate'); hold on;
% plot(er, B2, 'LineWidth', 1.5, 'DisplayName', 'Super strate'); 
xlim([4, 25]);
title('Relative BW vs. \epsilon_r');
ylabel('BW(%)');
xlabel('\epsilon_r');
grid on;
%legend show;

