%% This code takes in the necessary parameters for an exponential
%% fractional area distribution and outputs the scattering coefficient
%% under the assumption of independent single scattering

%% Version 2: This is an attempt to begin incorperating uncertainty into
%% the matlab code. Performing 10000 random samplings of rock
%% distributions, a standard deviation in the number of rocks per diameter
%% (binned into 1cm intervals) has been calculated in a seperate python
%% script. This will read in those values and see if the errors can
%% roughly match that of XFdtd.

%% Version 3: This is an attempt to including the covariance between the
%% number density variances. If you sample from a distribution until N
%% samples have been accepted, and you bin the different samples based on
%% their size, then there is a dependence on the number of samples in a 
%% particular bin based on the number of samples in all other bins.
%% Version 2 had assumed all the number densities were independent

function[kappa, Ekappa] = MieSphereIndependentScatV3(D, D_ave, k, alpha, beta, freq_GHz, eps_sp, eps_bk, DeltaD, WGSize)

q_k = alpha + beta/k; 

Es = zeros(length(D),length(freq_GHz)); Ea= Es; Eb=Es; Ee=Es; Qs = Es; 
n_vD = zeros(length(D), 1);

for iii = 1:length(freq_GHz)
for ii = 1:length(D)
    %% Ulaby's code calculates the efficiency factors using radii
    %% So have to remember to divide diameter by 2
    [Es(ii,iii), Ea(ii,iii), Ee(ii,iii), Eb(ii,iii)] = ...
        Mie_Rayleigh_ScatteringOfSpheres(D(ii)/2, freq_GHz(iii), conj(eps_sp), conj(eps_bk));
    Qs(ii, iii) = Es(ii, iii);

    %% We want the scattering cross section
    %% cross section = efficiency * pi * r^2
    Es(ii, iii) = Es(ii, iii)*pi*(D(ii)/2)^2;
    %n_vr(ii) = (k*q_k/(.8*pi*r_ave))*(exp(-2*q_k*rad(ii))/(rad(ii)^2));
    %r2n_vr(ii) = (2.8*k*q_k/(2*r_ave))*exp(-2*q_k*rad(ii));

    %% The Wu et al (CE4) rock distribution
    %% This is the volume number density of rocks with diameter D
    %% Paper reports Cumulative Fractional area. 
    %% See Golombek and Rapp (1997) for how to convert between cumulative frantional area
    %% and cumulative number of rocks with diameter greater than D
    %% Need volume density, so we extrapolate area to volume by dividing by an average diameter
    %% Peter's XF script used 0.11, so I've continued to use that but can easily be changed
    n_vD(ii) = (4*k*q_k/(pi*D_ave))*(exp(-q_k*D(ii))/(D(ii)^2));
end
end
Array = readtable(strcat('./std_n_VD_', WGSize, '10d.csv'), 'VariableNamingRule', 'preserve');
%Array = readtable("test_nVD_V2.csv", 'VariableNamingRule', 'preserve');
x = Array{:, 'STD(n_D)'};
%disp('number density');
%disp(x);
%disp('cross section');
%disp(Es(:, 10));
x2 = x.^2;  %% The dot means square each element in the vector
y = sqrt(n_vD);
%disp(transpose(y));
%disp(transpose(n_vD));
Es2 = transpose(Es).*transpose(Es); %% Same as above, the dot means perform the operation element wise

cov_freq = zeros(length(freq_GHz), 1);
disp("here");
%disp(transpose(Es(:,10)));
array = readtable(strcat('../RockNums_', WGSize, '10d_0.03D_0.064142Dm.csv'));
array = array{:,:};
array(:,1) = [];
cov_n = cov(transpose(array)/(8*8*10*0.01), 1);
disp(cov_n);
for ii = 1:length(cov_freq)
    cov_Es = Es(:,ii)*transpose(Es(:, ii));
    %cov_x = x*transpose(x) - diag(diag(x*transpose(x)));
    %disp(size(cov_Es*cov_x));
    %cov_freq(ii) = trace(cov_Es*cov_x);
    cov_freq(ii) = trace(cov_Es*cov_n);
end
disp('covariant term');
disp(sqrt(cov_freq));
%disp('covariant number density');
%disp(cov_x);
%disp(cov_freq);
%disp(Es2*x2);
y = sqrt(n_vD);
%% Scattering Coefficient is the integral of scattering cross section times volume number density over radius
kappa = transpose(Es)*n_vD*DeltaD;
disp(kappa);
%kappa = cov_freq*DeltaD;
Ekappa = sqrt(cov_freq)*DeltaD;
%disp(Ekappa);
%Ekappa = sqrt(Es2*x2)*DeltaD;
disp(sqrt(Es2*x2));
%Ekappa = transpose(Es)*x*DeltaD; %% While this solution is nonsensical (it
%results in a signal gain in attenuation) it does set an absolute upper and
%lower bound.

end