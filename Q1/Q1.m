clc
close all
clear all

load('Q1.mat');
fs = 100;
T = 100*fs;
t = 0.0001:1/fs:100;

figure('Name',"Main Signal");
for i = 1:8
    subplot(8,1,i);
    plot(t,X_org(i,:));
    if i==1 title("Main Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end

%% Part 1 GEVD
clc
%%%%% GEVD for periodic signal

% Main signal covariance
X = X_org - mean(X_org,2);
C_X = cov(X.');

% Periodic covariance
thau = 400;
X_shifted = circshift(X,thau,2);

for i = 1:8
    for j = 1:8
        P_X(i,j) = sum(X(i,:).*X_shifted(j,:)) / T;
    end
end
P_X = 0.5*(P_X + P_X.');

% GEVD
[W,Lambda] = eig(P_X,C_X);
[~,perm] = sort(diag(Lambda),'descend'); 
W = W(:,perm);
Lambda = Lambda(perm,perm);

% Finding sources and new X
Sources_1 = transpose(W)*X_org;

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_1(i,:));
    if i==1 title("Estimated Sources (period)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

A = inv(W.');
X_new = A(:,1)*Sources_1(1,:);


figure('Name',"Remaked Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new(i,:));
    if i==1 title("Remaked Signal (period)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X1(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new,X1));

%% Part 1 DSS
clc
% Whitening
X = X_org - mean(X_org,2);
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

figure('Name',"Whitened Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X(i,:));
    if i==1 title("Main Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8

    subplot(8,2,2*i);
    plot(t,Z(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

% Finding sources
W_1 = randn(8,1); % first W is normalized random
rrmse_old = 1000; % rrmse of former cycle
rrmse_new = 800; % rrmse of the new cycle
n = 0; % number of steps
while (rrmse_old-rrmse_new >0.001)
    n = n+1; % one more step
    r_1 = transpose(W_1)*Z; % estimating first source
    W = W_1; % holding W which was used in estimation

    % denoised source
    thau = 400;
    period_n = int32(T/thau);
    r = r_1(1+thau:2*thau);
    for i = 1:period_n
        r_1((i-1)*thau+1:i*thau) =r;
    end

    for j = 1:8
        W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
    end
    W_1 = W_1 / norm(W_1); % normalizing W
    Z_new = W*r_1; % finding new Z
    rrmse_old = rrmse_new;
    rrmse_new = RRMSE(Z_new,X1);
end
display("number of steps is : " + n);

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X1(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new,X2));


%% Part 2 GEVD
clc
%%%%% GEVD for unknown periodic signal
% Finding period
MAX = 0;
for thau = 300:700
    X_shifted = circshift(X,thau,2);
    for i = 1:8
        for j = 1:8
            c(i,j) = sum(X(i,:).*X_shifted(j,:)) / T;
        end
    end
    if (trace(c)>MAX)
        THAU = thau;
        MAX = trace(c);
    end
    c = zeros(8,8);
end
display("Best thau is : " + THAU);

% Periodic covariance
X_shifted = circshift(X,THAU,2);
for i = 1:8
    for j = 1:8
        P_X(i,j) = sum(X(i,:).*X_shifted(j,:)) / T;
    end
end
P_X = 0.5*(P_X + P_X.');

% GEVD
[W,Lambda] = eig(P_X,C_X);
[~,perm] = sort(diag(Lambda),'descend'); 
W = W(:,perm);
Lambda = Lambda(perm,perm);

% Finding sources and new X
Sources_1 = transpose(W)*X_org;

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_1(i,:));
    if i==1 title("Estimated Sources (UNK period)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

A = inv(W.');
X_new = A(:,1)*Sources_1(1,:);

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new(i,:));
    if i==1 title("Remade Signal (UNK period)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X1(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new,X1));


%% Part 2 DSS
clc
%%%%% DSS for unknown periodic signal
% finding best thau
X = X_org - mean(X_org,2);
MAX = 0;
for thau = 300:700
    X_shifted = circshift(X,thau,2);
    for i = 1:8
        for j = 1:8
            c(i,j) = sum(X(i,:).*X_shifted(j,:)) / T;
        end
    end
    if (trace(c)>MAX)
        THAU = thau;
        MAX = trace(c);
    end
    c = zeros(8,8);
end
display("Best thau is : " + THAU);

% Whitening
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

% Finding sources
W_1 = randn(8,1); % first W is normalized random
rrmse_old = 1000; % rrmse of former cycle
rrmse_new = 800; % rrmse of the new cycle
n = 0; % number of steps
while (rrmse_old-rrmse_new >0.001)
    n = n+1; % one more step
    r_1 = transpose(W_1)*Z; % estimating first source
    W = W_1; % holding W which was used in estimation

    % denoised source
    thau = 400;
    period_n = int32(T/thau);
    r = r_1(1+thau:2*thau);
    for i = 1:period_n
        r_1((i-1)*thau+1:i*thau) =r;
    end

    for j = 1:8
        W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
    end
    W_1 = W_1 / norm(W_1); % normalizing W
    Z_new = W*r_1; % finding new Z
    rrmse_old = rrmse_new;
    rrmse_new = RRMSE(Z_new,X1);
end
display("number of steps is : " + n);

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X1(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new,X2));


%% Part 3 GEVD
clc
%%%%% GEVD for burst-like signal

% Main signal covariance
X = X_org - mean(X_org,2);
C_X = cov(X.');

% Burst signal covariance
X_burst = T1.*X;
C_X_burst = cov(X_burst.');

% GEVD
[W,Lambda] = eig(C_X_burst,C_X);
[~,perm] = sort(diag(Lambda),'descend'); 
W = W(:,perm);
Lambda = Lambda(perm,perm);

% Finding sources and new X
Sources_2 = transpose(W)*X_org;

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_2(i,:));
    if i==1 title("Estimated Sources (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

A = inv(W.');
X_new = A(:,1)*Sources_2(1,:);

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new(i,:));
    if i==1 title("Remade Signal (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X2(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new,X2));


%% Part 3 DSS
clc
%%%%% DSS for for burst-like signal

% Whitening
X = X_org - mean(X_org,2);
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

% Finding sources
W_1 = randn(8,1); % first W is normalized random
rrmse_old = 100; % rrmse of former cycle
rrmse_new = 80; % rrmse of the new cycle
n = 0; % number of steps
while (rrmse_old-rrmse_new >0.001)
    n = n+1; % one more step
    r_1 = transpose(W_1)*Z; % estimating first source
    W = W_1; % holding W which was used in estimation
    r_1 = T1.*r_1; % denoised source
    for j = 1:8
        W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
    end
    W_1 = W_1 / norm(W_1); % normalizing W
    Z_new = W*r_1; % finding new Z
    rrmse_old = rrmse_new;
    rrmse_new = RRMSE(Z_new,X2);
end
display("number of steps is : " + n);

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X2(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new,X2));

%% Part 4 GEVD
clc
%%%%% GEVD for burst like signal

% Main signal covariance
X = X_org - mean(X_org,2);
C_X = cov(X.');

%%%% finding best of T2
MIN = 100;
for i = 1:1000
    T3 = ones(1,T) / i;
    T3 = T3 + T2;
    % Burst signal covariance
    X_burst_2 = T3.*X;
    C_X_burst_2 = cov(X_burst_2.');
    
    % GEVD
    [W_2,Lambda_2] = eig(C_X_burst_2,C_X);
    [~,perm] = sort(diag(Lambda_2),'descend'); 
    W_2 = W_2(:,perm);
    Lambda_2 = Lambda_2(perm,perm);
    
    % Finding sources and new X
    Sources_3 = transpose(W_2)*X_org;

    % Remaking signal
    A_2 = inv(W_2.');
    X_new_2 = A_2(:,1)*Sources_3(1,:);

    if (RRMSE(X_new_2,X2) < MIN)
        MIN = RRMSE(X_new_2,X2);
        X_new_best = X_new_2;
        Sources_best = Sources_3;
    end
end

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_best(i,:));
    if i==1 title("Estimated Sources (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new_best(i,:));
    if i==1 title("Remade Signal (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X2(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new_best,X2));


%% Part 4 DSS
%%%%% DSS for burst-like signal
clc

% Whitening
X = X_org - mean(X_org,2);
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

% Finding sources
W_1 = randn(8,1); % first W is normalized random
rrmse_old = 100; % rrmse of former cycle
rrmse_new = 80; % rrmse of the new cycle
n = 0; % number of steps
while (rrmse_old-rrmse_new >0.001)
    n = n+1; % one more step
    r_1 = transpose(W_1)*Z; % estimating first source
    W = W_1; % holding W which was used in estimation
    r_1 = T2.*r_1; % denoised source
    for j = 1:8
        W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
    end
    W_1 = W_1 / norm(W_1); % normalizing W
    Z_new = W*r_1; % finding new Z
    rrmse_old = rrmse_new;
    rrmse_new = RRMSE(Z_new,X2);
end
display("number of steps is : " + n);

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X2(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new,X2));


%% Part 5 GEVD
clc
%%%%% GEVD for bandpass signal

% Main signal covariance
X = X_org - mean(X_org,2);
C_X = cov(X.');

% Fourier transform
for i = 1:8
    X_f(i,:) = fftshift(fft(X(i,:)));
end
C_X_f = cov(X_f.');
C_X_f = 0.5*(C_X_f + C_X_f.');

% Filtering signal for 10-15 Hz
T3 = zeros(8,T);
for i = 6000:7500
    T3(:,i) = 1;
end
for i = 3500:4000
    T3(:,i) = 1;
end
X_f_filtered = X_f.*T3;
C_X_f_filtered = cov(X_f_filtered.');
C_X_f_filtered = 0.5*(C_X_f_filtered + C_X_f_filtered.');

% GEVD
[W,Lambda] = eig(C_X_f_filtered,C_X_f);
[~,perm] = sort(diag(Lambda),'descend'); 
W = W(:,perm);
Lambda = Lambda(perm,perm);

% Finding sources and new X
Sources_4 = transpose(W)*X_org;
A = inv(W.');
X_new = A(:,1)*Sources_4(1,:);

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_4(i,:));
    if i==1 title("Estimated Sources (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new(i,:));
    if i==1 title("Remade Signal (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X3(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new,X2));


%% Part 5 DSS
%%%%% DSS for bandpass signal
clc

% Whitening
X = X_org - mean(X_org,2);
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

% Finding sources
W_1 = randn(8,1); % first W is normalized random
rrmse_old = 100; % rrmse of former cycle
rrmse_new = 80; % rrmse of the new cycle
n = 0; % number of steps
while (rrmse_old-rrmse_new >0.001)
    n = n+1; % one more step
    r_1 = transpose(W_1)*Z; % estimating first source
    W = W_1; % holding W which was used in estimation

    %%% denoised source
    % Filtering signal for 10-15 Hz
    fpass = [10 15]; % Pass band
    [b,a] = butter(6,fpass/(fs/2), 'bandpass'); % butter filter
    r_1 = filtfilt(b,a,r_1); % applying filter

    for j = 1:8
        W_1(j,1) = sum(Z(j,:)*transpose(r_1)); % finding new W
    end
    W_1 = W_1 / norm(W_1); % normalizing W
    Z_new = W*r_1; % finding new Z
    rrmse_old = rrmse_new;
    rrmse_new = RRMSE(Z_new,X2);
end
display("number of steps is : " + n);

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X3(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new,X2));

%% Part 6 GEVD
clc
%%%%% GEVD for bandpass signal

% Main signal covariance
X = X_org - mean(X_org,2);
C_X = cov(X.');

% Fourier transform
for i = 1:8
    X_f(i,:) = fftshift(fft(X(i,:)));
end
C_X_f = cov(X_f.');
C_X_f = 0.5*(C_X_f + C_X_f.');

% Finding best bandpass in 5-25 Hz
L = 4;
MIN = 10000;
Bandpass = [0,0];
for i = 1:L
    for j = 50:(250-i*50)
        T3 = zeros(8,T);
        for x = (5000 + j*10):(5000 + j*10 + i*500)
            T3(:,x) = 1;
        end
        for x = (5000 - j*10 - i*500):(5000 - j*10)
            T3(:,x) = 1;
        end
        X_f_filtered = X_f.*T3;
        C_X_f_filtered = cov(X_f_filtered.');
        C_X_f_filtered = 0.5*(C_X_f_filtered + C_X_f_filtered.');
        
        % GEVD
        [W,Lambda] = eig(C_X_f_filtered,C_X_f);
        [~,perm] = sort(diag(Lambda),'descend'); 
        W = W(:,perm);
        Lambda = Lambda(perm,perm);
        
        % Finding sources and new X
        Sources_4 = transpose(W)*X_org;
        A = inv(W.');
        X_new = A(:,1)*Sources_4(1,:);
        %RRMSE(X_new,X3)

        if (RRMSE(X_new,X3) < MIN)
            MIN = RRMSE(X_new,X3);
            X_new_best = X_new;
            Sources_best = Sources_4;
            Bandpass = [j/10, j/10 + i*5];
        end
    end
end
display("Best band is : " + Bandpass(1) + " to " + Bandpass(2) + " Hz");

figure('Name',"Estimated Sources");
for i = 1:8
    subplot(8,1,i);
    plot(t,Sources_4(i,:));
    if i==1 title("Estimated Sources (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source number" + i);
    grid minor
end

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,X_new_best(i,:));
    if i==1 title("Remade Signal (burst)"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X3(i,:));
    if i==1 title("Ideal Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(X_new_best,X3));


%% Part 6 DSS
%%%%% DSS for burst-like signal
clc

% Whitening
X = X_org - mean(X_org,2);
cov_signal = cov(X.'); % Covariance matrix of X
[V,Lambda] = eig(cov_signal); % D is eigen value matrix and V is eigen vectors
D = diag(diag(Lambda).^(-0.5))* V.';
Z = D*X;

% Finding sources
MIN = 10000;
for i = 1:4
    for j = 5:25-i*5
        W_1 = randn(8,1); % first W is normalized random
        rrmse_old = 100; % rrmse of former cycle
        rrmse_new = 80; % rrmse of the new cycle
        n = 0; % number of steps
        while (rrmse_old-rrmse_new >0.001)
            n = n+1; % one more step
            r_1 = transpose(W_1)*Z; % estimating first source
            W = W_1; % holding W which was used in estimation
        
            %%% denoised source
            % Filtering signal
            fpass = [j j+(i*5)]; % Pass band
            [b,a] = butter(6,fpass/(fs/2), 'bandpass'); % butter filter
            r_1 = filtfilt(b,a,r_1); % applying filter
        
            for x = 1:8
                W_1(x,1) = sum(Z(x,:)*transpose(r_1)); % finding new W
            end
            W_1 = W_1 / norm(W_1); % normalizing W
            Z_new = W*r_1; % finding new Z
            rrmse_old = rrmse_new;
            rrmse_new = RRMSE(Z_new,X3);
        end
        if (RRMSE(Z_new,X3) < MIN)
            MIN = RRMSE(Z_new,X3);
            Z_new_best = Z_new;
            r_1_best = r_1;
            Bandpass = [j,j+(i*5)];
        end
    end
end
display("number of steps is : " + n);
display("Best band is : " + Bandpass(1) + " to " + Bandpass(2) + " Hz");

% Remaking signal
figure('Name',"Source found")
plot(t,r_1);
title("Source found");
xlabel("Time(s)");
grid minor

figure('Name',"Remade Signal");
for i = 1:8
    subplot(8,2,2*i-1);
    plot(t,Z_new_best(i,:));
    if i==1 title("Remade Signal"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Channel" + i);
    grid minor
end
for i = 1:8
    subplot(8,2,2*i);
    plot(t,X3(i,:));
    if i==1 title("Whitened Source"); end
    if i==8 xlabel("Time(s)"); end
    ylabel("Source" + i);
    grid minor
end

%RRMSE
display("RRMSE of new Signal is : " + RRMSE(Z_new_best,X3));

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