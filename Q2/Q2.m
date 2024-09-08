clc
close all
clear all

%% Part 1
clc
load('Ex2.mat');
load("Electrodes.mat");
fs = 200;
t = 0.0001:1/fs:51.2;
T = length(X_org(1,:));

% Plotting signal
offset = max(max(abs(X_org)))/3 ;
disp_eeg(X_org,offset,fs,Electrodes.labels);
title("Unnoisy EEG")
xlim('tight');
grid minor

figure('Name',"Channel T3")
plot(t,X_org);
title("Finding spikes' times")
xlabel("Time(s)");
xlim('tight');
grid minor;

% Defining T1 as an on off period
spiky_max_points = [2457 2617 2996 3942 4663 5515 5730 6437 8238 9319 10190];
Upper_thresh = 120;
Lower_thresh = 180;
T1 = zeros(1,T); % on or off vector
for i = 1:11
    if (i~=11)
        for j = spiky_max_points(i)-Lower_thresh:spiky_max_points(i)+Upper_thresh
            T1(1,j) = 1;
        end
    else
        for j = spiky_max_points(i)-Lower_thresh:T
            T1(1,j) = 1;
        end
    end
end

%% Part 2 Making noisy signals
%%%%% Calculating energies
p_sig = 0;
p_noise1 = 0;
p_noise2 = 0;
for i = 1:length(X_org(:,1))
    for j = 1:length(X_org(1,:))
        p_sig = p_sig + X_org(i,j)^2; % energy of main signal
        p_noise1 = p_noise1 + X_noise_1(i,j)^2; % energy of noise
        p_noise2 = p_noise2 + X_noise_3(i,j)^2; % energy of noise
    end
end

%%%%% Calculating Sigmas
snr = -10;
sigma_snr10_1 = ((p_sig/p_noise1)*10^(snr/-10))^0.5;
sigma_snr10_2 = ((p_sig/p_noise2)*10^(snr/-10))^0.5;
snr = -20;
sigma_snr20_1 = ((p_sig/p_noise1)*10^(snr/-10))^0.5;
sigma_snr20_2 = ((p_sig/p_noise2)*10^(snr/-10))^0.5;

%%%%% Producing noisy signals with given SNRs
X_10_1 = X_org + sigma_snr10_1*X_noise_1; % noisy signal with -10 SNR
X_20_1 = X_org + sigma_snr20_1*X_noise_1; % noisy signal with -20 SNR
X_10_2 = X_org + sigma_snr10_2*X_noise_3; % noisy signal with -10 SNR
X_20_2 = X_org + sigma_snr20_2*X_noise_3; % noisy signal with -20 SNR

%% Part 2 GEVD
clc
%%%%% GEVD for burst-like signals

% Noisy signals' covariances
X_10_1 = X_10_1 - mean(X_10_1,2);
X_20_1 = X_20_1 - mean(X_20_1,2);
X_10_2 = X_10_2 - mean(X_10_2,2);
X_20_2 = X_20_2 - mean(X_20_2,2);
C_X_10_1 = cov(X_10_1.');
C_X_20_1 = cov(X_20_1.');
C_X_10_2 = cov(X_10_2.');
C_X_20_2 = cov(X_20_2.');

% Burst signal covariance
X_10_1_burst = T1.*X_10_1;
X_20_1_burst = T1.*X_20_1;
X_10_2_burst = T1.*X_10_2;
X_20_2_burst = T1.*X_20_2;
C_X_10_1_burst = cov(X_10_1_burst.');
C_X_20_1_burst = cov(X_20_1_burst.');
C_X_10_2_burst = cov(X_10_2_burst.');
C_X_20_2_burst = cov(X_20_2_burst.');

% GEVD
[W_10_1,Lambda_10_1] = eig(C_X_10_1_burst,C_X_10_1);
[W_20_1,Lambda_20_1] = eig(C_X_20_1_burst,C_X_20_1);
[W_10_2,Lambda_10_2] = eig(C_X_10_2_burst,C_X_10_2);
[W_20_2,Lambda_20_2] = eig(C_X_20_2_burst,C_X_20_2);
[~,perm_10_1] = sort(diag(Lambda_10_1),'descend'); 
[~,perm_20_1] = sort(diag(Lambda_20_1),'descend'); 
[~,perm_10_2] = sort(diag(Lambda_10_2),'descend'); 
[~,perm_20_2] = sort(diag(Lambda_20_2),'descend'); 
W_10_1 = W_10_1(:,perm_10_1);
W_20_1 = W_20_1(:,perm_20_1);
W_10_2 = W_10_2(:,perm_10_2);
W_20_2 = W_20_2(:,perm_20_2);
Lambda_10_1 = Lambda_10_1(perm_10_1,perm_10_1);
Lambda_20_1 = Lambda_20_1(perm_20_1,perm_20_1);
Lambda_10_2 = Lambda_10_2(perm_10_2,perm_10_2);
Lambda_20_2 = Lambda_20_2(perm_20_2,perm_20_2);

% Finding sources and Spiky sources
Sources_10_1 = transpose(W_10_1)*X_10_1;
Sources_20_1 = transpose(W_20_1)*X_20_1;
Sources_10_2 = transpose(W_10_2)*X_10_2;
Sources_20_2 = transpose(W_20_2)*X_20_2;

SelSources_1 = [1 3 4];
SelSources_2 = [1 3 4 7];

% Denoised signals
A_10_1 = inv(W_10_1.');
A_20_1 = inv(W_20_1.');
A_10_2 = inv(W_10_2.');
A_20_2 = inv(W_20_2.');
X_new_10_1 = A_10_1(:,SelSources_1)*Sources_10_1(SelSources_1,:);
X_new_20_1 = A_20_1(:,SelSources_1)*Sources_20_1(SelSources_1,:);
X_new_10_2 = A_10_2(:,SelSources_2)*Sources_10_2(SelSources_2,:);
X_new_20_2 = A_20_2(:,SelSources_2)*Sources_20_2(SelSources_2,:);

% Plotting Sources
offset = max(max(abs(Sources_10_1)))/3 ;
disp_eeg(Sources_10_1,offset,fs);
title("Sources noise1 SNR -10")
xlim('tight');
grid minor
offset = max(max(abs(Sources_20_1)))/3 ;
disp_eeg(Sources_20_1,offset,fs);
title("Sources noise1 SNR -20")
xlim('tight');
grid minor
offset = max(max(abs(Sources_10_2)))/3 ;
disp_eeg(Sources_10_2,offset,fs);
title("Sources noise3 SNR -10")
xlim('tight');
grid minor
offset = max(max(abs(Sources_20_2)))/3 ;
disp_eeg(Sources_20_2,offset,fs);
title("Sources noise3 SNR -20")
xlim('tight');
grid minor

% Plotting Denoised Signals
offset = max(max(abs(X_new_10_1)))/3 ;
disp_eeg(X_new_10_1,offset,fs,Electrodes.labels);
title("Denoised EEG noise1 SNR -10")
xlim('tight');
grid minor
offset = max(max(abs(X_new_20_1)))/3 ;
disp_eeg(X_new_20_1,offset,fs,Electrodes.labels);
title("Denoised EEG noise1 SNR -20")
xlim('tight');
grid minor
offset = max(max(abs(X_new_10_2)))/3 ;
disp_eeg(X_new_10_2,offset,fs,Electrodes.labels);
title("Denoised EEG noise3 SNR -10")
xlim('tight');
grid minor
offset = max(max(abs(X_new_20_2)))/3 ;
disp_eeg(X_new_20_2,offset,fs,Electrodes.labels);
title("Denoised EEG noise3 SNR -20")
xlim('tight');
grid minor


%RRMSE
display("RRMSE of new Signal (noise1 SNR -10) is : " + RRMSE(X_new_10_1,X_org));
display("RRMSE of new Signal (noise1 SNR -20) is : " + RRMSE(X_new_20_1,X_org));
display("RRMSE of new Signal (noise3 SNR -10) is : " + RRMSE(X_new_10_2,X_org));
display("RRMSE of new Signal (noise3 SNR -20) is : " + RRMSE(X_new_20_2,X_org));

%% Part 2 DSS
clc
%%%%% DSS for burst-like signals

% Whitening
X_10_1 = X_10_1 - mean(X_10_1,2);
X_20_1 = X_20_1 - mean(X_20_1,2);
X_10_2 = X_10_2 - mean(X_10_2,2);
X_20_2 = X_20_2 - mean(X_20_2,2);
C_X_10_1 = cov(X_10_1.');
C_X_20_1 = cov(X_20_1.');
C_X_10_2 = cov(X_10_2.');
C_X_20_2 = cov(X_20_2.');
[V_10_1,Lambda_10_1] = eig(C_X_10_1); % D is eigen value matrix and V is eigen vectors
[V_20_1,Lambda_20_1] = eig(C_X_20_1);
[V_10_2,Lambda_10_2] = eig(C_X_10_2);
[V_20_2,Lambda_20_2] = eig(C_X_20_2);
D_10_1 = diag(diag(Lambda_10_1).^(-0.5))* V_10_1.';
D_20_1 = diag(diag(Lambda_20_1).^(-0.5))* V_20_1.';
D_10_2 = diag(diag(Lambda_10_2).^(-0.5))* V_10_2.';
D_20_2 = diag(diag(Lambda_20_2).^(-0.5))* V_20_2.';
Z_10_1 = D_10_1*X_10_1;
Z_20_1 = D_20_1*X_20_1;
Z_10_2 = D_10_2*X_10_2;
Z_20_2 = D_20_2*X_20_2;

%%%% Finding sources

%%% Noise1 SNR -10
Z = Z_10_1;
N = 3; % Number of sources
r = zeros(N,T);
steps = zeros;
for i = 1:N
    W_1 = randn(32,1);
    rrmse_old = 100; % rrmse of former cycle
    rrmse_new = 80; % rrmse of the new cycle
    n = 0; % number of steps
    while (rrmse_old-rrmse_new >0.0001)
        n = n+1; % one more step
        r_1 = transpose(W_1)*Z; % estimating first source
        W = W_1; % holding W which was used in estimation
        r_1 = T1.*r_1; % denoised source
        for j = 1:32
            W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
        end
        W_1 = W_1 / norm(W_1); % normalizing W
        Z_new_10_1 = W*r_1; % finding new Z
        rrmse_old = rrmse_new;
        rrmse_new = RRMSE(Z_new_10_1,X_org);
    end
    r(i,:) = r_1;
    Z = Z-Z_new_10_1;
    steps(i) = n;
end
%display("number of steps is : " + steps(:));
%display("RRMSE of new Signal is : " + RRMSE(Z_new_10_1,X_org));

% Plotting Sources
figure('Name',"Sources");
for i = 1:N
    subplot(N,1,i);
    plot(t,r(i,:));
    if i==1 title("Sources noise1 SNR -10"); end
    if i==N xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

% Plotting Denoised Signal
offset = max(max(abs(Z_new_10_1)))/3 ;
disp_eeg(Z_new_10_1,offset,fs,Electrodes.labels);
title("Denoised EEG noise1 SNR -10")
xlim('tight');
grid minor


%%% Noise1 SNR -20
Z = Z_20_1;
N = 3; % Number of sources
r = zeros(N,T);
steps = zeros;
for i = 1:N
    W_1 = randn(32,1);
    rrmse_old = 100; % rrmse of former cycle
    rrmse_new = 80; % rrmse of the new cycle
    n = 0; % number of steps
    while (rrmse_old-rrmse_new >0.0001)
        n = n+1; % one more step
        r_1 = transpose(W_1)*Z; % estimating first source
        W = W_1; % holding W which was used in estimation
        r_1 = T1.*r_1; % denoised source
        for j = 1:32
            W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
        end
        W_1 = W_1 / norm(W_1); % normalizing W
        Z_new_20_1 = W*r_1; % finding new Z
        rrmse_old = rrmse_new;
        rrmse_new = RRMSE(Z_new_20_1,X_org);
    end
    r(i,:) = r_1;
    Z = Z-Z_new_20_1;
    steps(i) = n;
end
%display("number of steps is : " + steps(:));
%display("RRMSE of new Signal is : " + RRMSE(Z_new_20_1,X_org));

% Plotting Sources
figure('Name',"Sources");
for i = 1:N
    subplot(N,1,i);
    plot(t,r(i,:));
    if i==1 title("Sources noise1 SNR -20"); end
    if i==N xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

% Plotting Denoised Signal
offset = max(max(abs(Z_new_20_1)))/3 ;
disp_eeg(Z_new_20_1,offset,fs,Electrodes.labels);
title("Denoised EEG noise1 SNR -20")
xlim('tight');
grid minor


%%% Noise3 SNR -10
Z = Z_10_2;
N = 2; % Number of sources
r = zeros(N,T);
steps = zeros;
for i = 1:N
    W_1 = randn(32,1);
    rrmse_old = 100; % rrmse of former cycle
    rrmse_new = 80; % rrmse of the new cycle
    n = 0; % number of steps
    while (rrmse_old-rrmse_new >0.0001)
        n = n+1; % one more step
        r_1 = transpose(W_1)*Z; % estimating first source
        W = W_1; % holding W which was used in estimation
        r_1 = T1.*r_1; % denoised source
        for j = 1:32
            W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
        end
        W_1 = W_1 / norm(W_1); % normalizing W
        Z_new_10_2 = W*r_1; % finding new Z
        rrmse_old = rrmse_new;
        rrmse_new = RRMSE(Z_new_10_2,X_org);
    end
    r(i,:) = r_1;
    Z = Z-Z_new_10_2;
    steps(i) = n;
end
%display("number of steps is : " + steps(:));
display("RRMSE of new Signal is : " + RRMSE(Z_new_10_2,X_org));

% Plotting Sources
figure('Name',"Sources");
for i = 1:N
    subplot(N,1,i);
    plot(t,r(i,:));
    if i==1 title("Sources noise3 SNR -10"); end
    if i==N xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

% Plotting Denoised Signal
offset = max(max(abs(Z_new_10_2)))/3 ;
disp_eeg(Z_new_10_2,offset,fs,Electrodes.labels);
title("Denoised EEG noise3 SNR -10")
xlim('tight');
grid minor


%%% Noise1 SNR -20
Z = Z_20_2;
N = 2; % Number of sources
r = zeros(N,T);
steps = zeros;
for i = 1:N
    W_1 = randn(32,1);
    rrmse_old = 100; % rrmse of former cycle
    rrmse_new = 80; % rrmse of the new cycle
    n = 0; % number of steps
    while (rrmse_old-rrmse_new >0.0001)
        n = n+1; % one more step
        r_1 = transpose(W_1)*Z; % estimating first source
        W = W_1; % holding W which was used in estimation
        r_1 = T1.*r_1; % denoised source
        for j = 1:32
            W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
        end
        W_1 = W_1 / norm(W_1); % normalizing W
        Z_new_20_2 = W*r_1; % finding new Z
        rrmse_old = rrmse_new;
        rrmse_new = RRMSE(Z_new_20_2,X_org);
    end
    r(i,:) = r_1;
    Z = Z-Z_new_20_2;
    steps(i) = n;
end
%display("number of steps is : " + steps(:));
%display("RRMSE of new Signal is : " + RRMSE(Z_new_20_2,X_org));

% Plotting Sources
figure('Name',"Sources");
for i = 1:N
    subplot(N,1,i);
    plot(t,r(i,:));
    if i==1 title("Sources noise3 SNR -20"); end
    if i==N xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

% Plotting Denoised Signal
offset = max(max(abs(Z_new_20_2)))/3 ;
disp_eeg(Z_new_20_2,offset,fs,Electrodes.labels);
title("Denoised EEG noise3 SNR -20")
xlim('tight');
grid minor


%% Part 3
% denoised1 -10 SNR
figure("Name","Part3-10-1");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-1 (-10)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,X_new_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel GEVD x-den-1')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,Z_new_10_1(13,:));
grid minor;
xlim('tight');
title('13th channel DSS x-den-1')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-1 (-10)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,X_new_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel GEVD x-den-1')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,Z_new_10_1(24,:));
grid minor;
xlim('tight');
title('24th channel DSS x-den-1')
xlabel('Time (s)')

% denoised1 -20 SNR
figure("Name","Part3-20-1");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-1 (-20)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,X_new_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel GEVD x-den-1')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,Z_new_20_1(13,:));
grid minor;
xlim('tight');
title('13th channel DSS x-den-1')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-1 (-20)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,X_new_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel GEVD x-den-1')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,Z_new_20_1(24,:));
grid minor;
xlim('tight');
title('24th channel DSS x-den-1')
xlabel('Time (s)')

% denoised2 -10 SNR
figure("Name","Part3-10-2");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-2 (-10)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,X_new_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel GEVD x-den-2')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,Z_new_10_2(13,:));
grid minor;
xlim('tight');
title('13th channel DSS x-den-2')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-2 (-10)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,X_new_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel GEVD x-den-2')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,Z_new_10_2(24,:));
grid minor;
xlim('tight');
title('24th channel DSS x-den-2')
xlabel('Time (s)')

% denoised2 -20 SNR
figure("Name","Part3-20-2");
subplot(4,2,1)
plot(t,X_org(13,:));
grid minor;
xlim('tight');
title('13th channel x-org')
xlabel('Time (s)')

subplot(4,2,3)
plot(t,X_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel x-noise-2 (-20)')
xlabel('Time (s)')

subplot(4,2,5)
plot(t,X_new_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel GEVD x-den-2')
xlabel('Time (s)')

subplot(4,2,7)
plot(t,Z_new_20_2(13,:));
grid minor;
xlim('tight');
title('13th channel DSS x-den-2')
xlabel('Time (s)')

subplot(4,2,2)
plot(t,X_org(24,:));
grid minor;
xlim('tight');
title('24th channel x-org')
xlabel('Time (s)')

subplot(4,2,4)
plot(t,X_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel x-noise-2 (-20)')
xlabel('Time (s)')

subplot(4,2,6)
plot(t,X_new_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel GEVD x-den-2')
xlabel('Time (s)')

subplot(4,2,8)
plot(t,Z_new_20_2(24,:));
grid minor;
xlim('tight');
title('24th channel DSS x-den-2')
xlabel('Time (s)')

%% Part 4
clc
display("RRMSE of new Signal (GEVD) (noise1 & SNR -10) is : " + RRMSE(X_new_10_1,X_org));
display("RRMSE of new Signal (GEVD) (noise1 & SNR -20) is : " + RRMSE(X_new_20_1,X_org));
display("RRMSE of new Signal (GEVD) (noise3 & SNR -10) is : " + RRMSE(X_new_10_2,X_org));
display("RRMSE of new Signal (GEVD) (noise3 & SNR -20) is : " + RRMSE(X_new_20_2,X_org));
display("RRMSE of new Signal (DSS) (noise1 & SNR -10) is : " + RRMSE(Z_new_10_1,X_org));
display("RRMSE of new Signal (DSS) (noise1 & SNR -20) is : " + RRMSE(Z_new_20_1,X_org));
display("RRMSE of new Signal (DSS) (noise3 & SNR -10) is : " + RRMSE(Z_new_10_2,X_org));
display("RRMSE of new Signal (DSS) (noise3 & SNR -20) is : " + RRMSE(Z_new_20_2,X_org));

%% FUNCTIONS
function c = RRMSE(X,X_org)
    b = 0;
    a = 0;
    for i = 1:length(X_org(:,1))
        for j = 1:length(X_org(1,:))
            a = a + (X_org(i,j)-X(i,j))^2;
            b = b + (X_org(i,j))^2;
        end
    end
    c = (a/b)^0.5;
end

