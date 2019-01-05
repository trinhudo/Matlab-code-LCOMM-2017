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
pN = 1/5; % power allocation coefficient
pF = 4/5;
Rth_near = 1; % bits/s/Hz
g1 = 2^(2*Rth_near)-1; % gamma_1
Rth_far = 1;
g2 = 2^(2*Rth_far)-1; % gamma_2
g2_non = 2^Rth_far - 1;
SimTimes = 5*10^0; % Monte-Carlo repetitions
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
    disp(strcat('SNR=',num2str(snravg_dB(ss))));
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
    % find the best near and far ones
    [gSNs(:,1),Nbest] = max(gSNi,[],2);
    [gSFs(:,1),Fbest] = max(gSFj,[],2);
    for yy = 1:SimTimes
        gNsFs(yy,1,1) = gNiFj(yy,Nbest(yy),Fbest(yy));
    end
    % SNR modelling
    snrSNs_xFs = (1-rho).*pF.*snravg(ss).*gSNs./...
        ((1-rho).*pN.*snravg(ss).*gSNs + 1);
    %
    snrSNs_xNs = (1-rho).*pN.*snravg(ss).*gSNs;
    snrSFs = pF.*snravg(ss).*gSFs./...
        (pN.*snravg(ss).*gSFs + 1);
    snrNsFs_DF = rho.*eta.*snravg(ss).*gSNs.*gNsFs;
    snrNsFs_AF = rho.*(1-rho).*eta.*(snravg(ss).^2).*pF.*(gSNs.^2).*gNsFs./...
        (rho.*(1-rho).*eta.*(snravg(ss).^2).*pN.*(gSNs.^2).*gNsFs + ...
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
        %% for DF only
        if (snrSNs_xFs(zz) >= g2) && (max(snrSFs(zz),snrNsFs_DF(zz))<g2)
            countDF = countDF + 1;
        elseif (snrSNs_xFs(zz) < g2) && (snrSFs(zz) < g2)
            countDF = countDF + 1;
        end
        %% for AF only
        if (max(snrSFs(zz),snrNsFs_AF(zz)) < g2)
            countAF = countAF + 1;
        end
    end
    PoutHybrid(ss) = countHybrid/SimTimes;
    PoutDF(ss) = countDF/SimTimes;
    PoutAF(ss) = countAF/SimTimes;
    %% Analytical results
    a1 = (1-rho)*pF*snravg(ss);
    a2 = (1-rho)*pN*snravg(ss);
    b1 = pF*snravg(ss);
    b2 = pN*snravg(ss);
    c  = rho*eta*snravg(ss);
    mu = g2/(a1-a2*g2);
    %
    Phi1A = 0;
    %
    for hh = 0:N
        Phi1A_temp = nchoosek(N,hh)*((-1)^hh)*...
            exp(-hh*lSN*g2/(a1-a2*g2));
        Phi1A = Phi1A + Phi1A_temp;
    end
    %
    Phi1A_asym = (lSN*g2/(a1-a2*g2))^N;
    %
    Phi1B = 0;
    %
    for kk = 0:M
        Phi1B_temp = nchoosek(M,kk)*((-1)^kk)*...
            exp(-kk*lSF*g2/(b1-b2*g2));
        Phi1B = Phi1B + Phi1B_temp;
    end
    % asymptotic Phi1
    Phi1B_asym = (lSF*g2/(b1-b2*g2))^M;
    %%
    Phi2A = 0;
    Phi2A_asym_temp = 0;
    Phi2A_asym = 0;
    %
    for jj=1:N
        %
        Phi2A_temp = nchoosek(N,jj)*((-1)^(jj+1))*...
            jj*lSN*lNF*g2/c*(-ei(-jj*lSN*mu));
        Phi2A = Phi2A + Phi2A_temp;
        %
        Phi2_temp_asym = nchoosek(N,jj)*((-1)^(jj+1))*...
            jj*lSN*lNF*g2/rho/eta/snravg(ss)*...
            (-log(jj*lSN*mu));
        Phi2A_asym_temp = Phi2A_asym_temp + Phi2_temp_asym;
        Phi2A_asym = Phi2A_asym_temp;
    end
    %%
    Psi_sum = 0;
    Psi_sum_asym = 0;
    %
    for ll = 1:N
        tt = 4*ll*lSN*lNF*snravg(ss)*mu/c;
        %
        Psi_temp = nchoosek(N,ll)*((-1)^(ll+1))*...
            exp(-ll*lSN*mu)*...
            sqrt(tt)*besselk(1,sqrt(tt));
        Psi_sum = Psi_sum + Psi_temp;
        %
        Psi_sum_temp_asym = nchoosek(N,ll)*((-1)^(ll-1))*...
            tt/2*log(sqrt(tt)/2);
        Psi_sum_asym = Psi_sum_asym + Psi_sum_temp_asym;
    end
    %
    Psi = 1 - Psi_sum;
    Psi_asym = - Psi_sum_asym;
    %
    Theta = 0;
    for vv = 1:N
        Theta_temp = nchoosek(N,vv)*((-1)^(vv+1))*...
            (1-exp(-vv*lSN*g2/(a1-a2*g2)));
        Theta = Theta + Theta_temp;
    end
    %
    Theta_asym = (lSN*g2/(a1-a2*g2))^N;
    % Pout of DF
    PoutDFana(ss) = Phi1B*(Phi1A+Phi2A);
    PoutAFana(ss) = Phi1B*Psi;
    PoutHybridana(ss) = Phi1B*(Phi2A+Theta);
    %
    PoutDF_asym(ss) = Phi1B_asym*(Phi1A_asym+Phi2A_asym);
    PoutAF_asym(ss) = Phi1B_asym*Psi_asym;
    PoutHybrid_asym(ss) = Phi1B_asym*(Phi2A_asym+Theta_asym);
end
%% plot
Ns_sim = [0.00970510000000000,0.000383900000000000,1.29000000000000e-05,4.300000000000000e-07,0,0];
Ns_ana = [0.00971433625879054,0.000390396083015721,1.33361475804783e-05,4.32207955980424e-07,1.37743330075324e-08,4.36655822610987e-10];
Ns_asym = [0.0138240000000000,0.000437153263741677,1.38240000000000e-05,4.37153263741678e-07,1.38240000000000e-08,4.37153263741678e-10];
%
PoutHybrid = [0.0372135200000000,0.00954344000000000,0.00142006000000000,5.15400000000000e-05,7.40000000000000e-07,0];
PoutDF = [0.0372135200000000,0.00954344000000000,0.00142006000000000,5.15400000000000e-05,7.40000000000000e-07,0];
PoutAF = [0.300254820000000,0.0953094200000000,0.0143643400000000,0.000516880000000000,7.82000000000000e-06,1.150000000000000e-07];
%
hplot1= semilogy(snravg_dB,Ns_sim,'kd',...
    snravg_dB,PoutHybrid,'ro',...
    snravg_dB,PoutDF,'b+',...
    snravg_dB,PoutAF,'g>',...
    snravg_dB,Ns_ana,'k-',...
    snravg_dB,Ns_asym,'k:',...
    snravg_dB,PoutHybridana,'r-',...
    snravg_dB,PoutHybrid_asym,'r:');
hold on
%
hplot2 = semilogy(snravg_dB,PoutDFana,'b-',...
    snravg_dB,PoutDF_asym,'b:',...
    snravg_dB,PoutAFana,'g-',...
    snravg_dB,PoutAF_asym,'g:');
%
legend('N_s, (sim.)','F_s, Hybrid (sim.)','F_s, DF (sim.)','F_s, AF (sim.)')
%
xlabel('Transmit SNR (dB)')
ylabel('Outage Probability')
%
axis([0 25 10^-8 1])
set(gca,'XTick',0:5:25)
%
% legend(hplot1, 'N_s (sim.)','N_s (ana.)','N_s (asym.)',...
%     'F_s, Hybrid (sim.)','F_s, Hybrid (ana.)','F_s, Hybrid (asym.)', 'Location','SouthWest'); %display legend 1
% ax=axes('Position',get(gca,'Position'),'Visible','Off');
% legend(ax, hplot2, 'F_s, DF (sim.)','F_s, DF (ana.)','F_s, DF (asym.)',...
%     'F_s, AF (sim.)','F_s, AF (ana.)','F_s, AF (asym.)', 'Location','NorthEast');