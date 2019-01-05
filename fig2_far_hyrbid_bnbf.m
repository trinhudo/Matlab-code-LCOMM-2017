format long
close all
clear all
%
%% Network topology
%
S = [0,0.5];
N_user = [0.2,0.5]; % near users
F_user = [1,0.5]; % far users
%
dSN = sqrt(abs(S(1)-N_user(1))^2 + abs(S(2)-N_user(2))^2);
dSF = sqrt(abs(S(1)-F_user(1))^2 + abs(S(2)-F_user(2))^2);
dNF = sqrt(abs(F_user(1)-N_user(1))^2 + abs(F_user(2)-N_user(2))^2);
%
%% Simulation parameters
%
N = 3; % # of near users
M = 3; % # of far users
rho = 0.5; % power splitting ratio
snravg_dB = 0:5:25; % transmit SNR = Ps/N0 in dB
snravg = 10.^(snravg_dB./10);
epsilon = 3; % pathloss exponent
% Average channel gain: (lambda) = 1/(d^(-path_loss_exponent))
lSN = dSN^epsilon; % lambda
lSF = dSF^epsilon;
lNF = dNF^epsilon;
%
eta = 0.7; % energy conversion coefficient
thetaN = 1/5; % power allocation coefficient
thetaF = 4/5;
Rth_near = 1; % bits/s/Hz
g1 = 2^(2*Rth_near)-1; % gamma_1
Rth_far = 1;
g2 = 2^(2*Rth_far)-1; % gamma_2
g2_non = 2^Rth_far - 1;
SimTimes = 2*10^7; % Monte-Carlo repetitions
%
%% Simulation
%
hSFj = zeros(SimTimes,M);
hSNi = zeros(SimTimes,N);
hNiFj = zeros(SimTimes,N,M);
gNsFs = zeros(SimTimes,1,1);
PoutHybrid = zeros(1,length(snravg_dB));
%
for ss = 1:length(snravg_dB)
    disp(strcat('SNR=',...
        num2str(snravg_dB(ss))));
    % channel modelling
    for kk = 1: M
        hSFj(:,kk) = sqrt(1/2/lSF)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
    end
    for nn = 1:N
        hSNi(:,nn) = sqrt(1/2/lSN)*...
            (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        for kk = 1:M
            hNiFj(:,nn,kk) = sqrt(1/2/lNF)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
        end
    end
    % channel gains
    gSNi = abs(hSNi.^2);
    gSFj = abs(hSFj.^2);
    gNiFj = abs(hNiFj.^2);
    % Find the best near and far ones
    %     [gSNs(:,1),Nbest] = max(gSNi,[],2);
    %     [gSFs(:,1),Fbest] = min(gSFj,[],2);
    %     for yy = 1:SimTimes
    %         gNsFs(yy,1,1) = gNiFj(yy,Nbest(yy),Fbest(yy));
    %     end
    % Randomly selection
    for yy = 1:SimTimes
        Near_rand = randperm(N,1);
        Far_rand = randperm(M,1);
        %
        gSNs(yy,1) = gSNi(yy,Near_rand);
        gSFs(yy,1) = gSFj(yy,Far_rand);
        gNsFs(yy,1,1) = gNiFj(yy,Near_rand,Far_rand);
    end
    % SNR modelling
    snrSNs_xFs = (1-rho).*thetaF.*snravg(ss).*gSNs./...
        ((1-rho).*thetaN.*snravg(ss).*gSNs + 1);
    %
    snrSNs_xNs = (1-rho).*thetaN.*snravg(ss).*gSNs;
    snrSFs = thetaF.*snravg(ss).*gSFs./...
        (thetaN.*snravg(ss).*gSFs + 1);
    snrNsFs_DF = rho.*eta.*snravg(ss).*gSNs.*gNsFs;
    snrNsFs_AF = rho.*(1-rho).*eta.*(snravg(ss).^2).*thetaF.*(gSNs.^2).*gNsFs./...
        (rho.*(1-rho).*eta.*(snravg(ss).^2).*thetaN.*(gSNs.^2).*gNsFs + ...
        rho.*eta.*snravg(ss).*gSNs.*gNsFs + snravg(ss).*gSNs + 1);
    % count outage events
    countHybrid = 0;
    countDF = 0;
    countAF = 0;
    %
    for zz = 1:SimTimes
        %% Hybrid
        if (snrSNs_xFs(zz) >= g2) && (max(snrSFs(zz),snrNsFs_DF(zz)) < g2)
            countHybrid = countHybrid + 1;
        elseif (snrSNs_xFs(zz) < g2) && (max(snrSFs(zz),snrNsFs_AF(zz)) < g2)
            countHybrid = countHybrid + 1;
        end
    end
    PoutHybrid(ss) = countHybrid/SimTimes;
end
%% plot
semilogy(snravg_dB,PoutHybrid,'r+-')

% h=legend('Hybrid BNBF (sim.)');

xlabel('Transmit SNR (dB)')
ylabel('Outage Probability')
%
axis([0 25 10^-8 1])
set(gca,'XTick',0:5:25)
