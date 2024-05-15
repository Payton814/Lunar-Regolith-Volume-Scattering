clear;

freq_GHz = 0:0.001:0.9; % input frequency in GHz
%rad = 0.05:0.001:0.625; % input radius size in m
DeltaD = 0.01;
D = 0.03:DeltaD:1.24;
%disp(D);
WGSize = '4a4b'; %% This defines the waveguide size we want to look at


k = 0.0021;
alpha = 0.5648;
beta = 0.01258;
%q_k = 0.5648 + 0.01258/k;
%D_ave = 0.11;
%D_ave = 0.094519;
D_ave = 0.064142;

eps_bk = 3.01 + 1i * 0.0; % permitivity of the background material, note the imaginary part is 0 making the material lossless
eps_sp = 6.01 + 1i * 0.0; % permitivity of spherical scatterer

[kappa, Ekappa] = MieSphereIndependentScatV2(D, D_ave, k, alpha, beta, freq_GHz,eps_sp, eps_bk, DeltaD, WGSize);
%[kappa2] = MieSphereIndependentScat(D, D_ave, 0.00345, -10, 0.04558, freq_GHz,eps_sp, eps_bk);

Ppf_0 = ones(length(freq_GHz), 1); % initial Power per frequnecy
d = 10; % distance wave travels in meters
Ppf = zeros(length(freq_GHz), 1);
PpfUB = zeros(length(freq_GHz), 1);
PpfLB = zeros(length(freq_GHz), 1);
%Ppf2 = zeros(length(freq_GHz), 1);
%% The attenuation of a wave traveling through a medium with scattering coefficient ks
%% is exponential with distance, atten = exp(-ks*d)
for ii = 1:length(freq_GHz)
    Ppf(ii) = Ppf_0(ii)*exp(-kappa(ii)*d);
    PpfUB(ii) = Ppf_0(ii)*exp(-(kappa(ii)+Ekappa(ii))*d);
    PpfLB(ii) = Ppf_0(ii)*exp(-(kappa(ii)-Ekappa(ii))*d);
    %Ppf2(ii) = Ppf_0(ii)*exp(-kappa2(ii)*d);
end
disp(kappa - Ekappa)



files = dir(strcat('./CE4/', WGSize, '10d/*.csv')); %% Grabs all the CSV files in the specified folder and puts their names into an array called files
disp(size(files));
Trials = []; %% Initialize an empty array to be populated
for i=1:length(files)
    array = readtable(strcat('./CE4/', WGSize, '10d/', files(i).name), 'VariableNamingRule', 'preserve');
    x = array{150:3977, 'Frequency (GHz)'};
    y = array{150:3977, 'S2'};
    Trials = [Trials 10*log10(y(:))]; %% Concatinate the y value from the array that was just read in
end

subplot(1,1,1)
hold on
shadedErrorBar(x, transpose(Trials), {@mean,@std}, 'lineprops', '-b');

plot(freq_GHz, 10*log10(Ppf), '-k');
plot(freq_GHz, 10*log10(PpfUB), '-r');
plot(freq_GHz, 10*log10(PpfLB), '-r');
hold off
xlabel('Frequency (GHz)')
ylabel('Attenuation per frequency [dB]')
title('Attenuation per Frequency after traveling 10 m')
legend({'XFdtd 7x7x10', 'Matlab'})
%legend
grid on