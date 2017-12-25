function power_estimation_RPE
%POWER_ESTIMATION_RPE   Estimation of statistical power.
%   POWER_ESTIMATION_RPE estimates statistical power depending on effect
%   size and sample size to help determine the necessary sample size to
%   reach a desired level of statistical power.
%
%   POWER_ESTIMATION_RPE is developed to estimate the potential modulation
%   of firing rates of basal forebrain neurons by reward expectations or
%   reward prediction error. For estimating the parameters and the expected
%   effect size we used data from Hangya et al., 2015.
%
%   Briefly, we set 10%, 50% and 100% effect sizes as ‘small’, ‘medium’ and
%   ‘large’ effect. As reward often elicited a single spike per trial in
%   HDB cholinergic neurons with short latency (27 ms), small jitter (11
%   ms) and variable reliability depending on the signal-to-noise ratio of
%   the reward predicting cue (average, 0.45), we modeled the expected
%   effect using these parameters on top of a 5 Hz baseline Poisson firing
%   (Fig.1G in Hangya et al., 2015). We estimated that 60% of HDB neurons
%   may show a modulation and assumed no modulation in the remaining 40%;
%   nevertheless this estimation probably represents a lower bound, since
%   modulation of firing rate by surprise or reward expectation may not
%   reach statistical significance in all recorded neurons. To be on the
%   conservative side, we modeled a relatively short session with only 100
%   rewarded trials. Statistical significance was investigated by
%   performing an ROC analysis in a 50 ms window and tested with Wilcoxon
%   sign-rank test at alpha = 0.05. A hundred experiments were simulated to
%   provide an estimate of statistical power, that is, the proportion of
%   simulated experiments with statistically significant effect.
%
%   A description of the results can be found at
%   www.github.com/hangyabalazs.
%
%   Reference: Hangya B, Ranade SP, Lorenc M, Kepecs A (2015) Central
%   cholinergic neurons are rapidly recruited by reinforcement feedback.
%   Cell, 162:1155–1168.

%   Balazs Hangya
%   Laboratory of Systems Neuroscience, 
%   Hungarian Academy os Sceinces
%   hangya.balazs@koki.mta.hu

% Power estimation parameters
effect_sizes = [0.1 0.5 1];  % effect size
effect_probability = 0.6;   % not all cells may show the effect
sample_sizes = 5:20;   % sample size

% Power
for effect_size = effect_sizes
    next = 0;
    P = [];
    for sample_size = sample_sizes
        disp(sample_size)
        next = next + 1;
        [P(next) p H] = main(effect_size,effect_probability,sample_size); %#ok<AGROW>
    end
    switch find(ismember(effect_sizes,effect_size))
        case 1
            P_small = P;
        case 2
            P_medium = P;
        case 3
            P_large = P;
    end
end

% Plot
figure
plot(sample_sizes,P_small,'Color',[0.1 0.9 0.1],'LineWidth',2)
hold on
plot(sample_sizes,P_medium,'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(sample_sizes,P_large,'Color',[0.9 0.1 0.9],'LineWidth',2)
legend({'small' 'medium' 'large'})
xlabel('Sample size')
ylabel('Statistical power')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Power estimation parameters
effect_probability = 0.6;   % not all cells may show the effect
sample_sizes = [10 15 20];   % number of cells
effect_sizes = 0.1:0.1:1;   % effect size

% Power
for sample_size = sample_sizes
    next = 0;
    P = [];
    for effect_size = effect_sizes
        disp(effect_size)
        next = next + 1;
        [P(next) p H] = main(effect_size,effect_probability,sample_size); %#ok<AGROW>
    end
    switch find(ismember(sample_sizes,sample_size))
        case 1
            P_small = P;
        case 2
            P_medium = P;
        case 3
            P_large = P;
    end
end

% Plot
figure
plot(effect_sizes,P_small,'Color',[0.1 0.9 0.1],'LineWidth',2)
hold on
plot(effect_sizes,P_medium,'Color',[0.5 0.5 0.5],'LineWidth',2)
plot(effect_sizes,P_large,'Color',[0.9 0.1 0.9],'LineWidth',2)
legend({'small' 'medium' 'large'})
xlabel('Effect size')
ylabel('Statistical power')
axis square
keyboard

% -------------------------------------------------------------------------
function [P p H] = main(effect_size,effect_probability,sample_size)

% Simulation parameters
mx = 1000000;   % data length in ms (30 min)
esp = 10000;    % event spacing in ms
baselinerate = 5;   % baseline spike rate (Poisson) in Hz
esp_lat = 27; % latency of evoked spikes in ms (Gaussian)
esp_jit = 11; % jitter of evoked spikes in ms (Gaussian)
nesp = 1;     % number of evoked spikes per event (can be randomized later)
pesp = 0.45;   % response probability
wn = [-200 600];   % window (ms)
sno = 100;   % number of simulations

% ROC for simulated effect
dsp = false;
[p H] = deal(nan(1,sno));
for iS = 1:sno
%     disp(iS)
    [ROC P SE] = deal(nan(1,sample_size));
    inx1 = - wn(1);
    inx2 = inx1 + 50;   % window for ROC: 0-50 ms
    for iC = 1:sample_size
%         disp(iC)
        spt0 = spike_trains(mx,esp,baselinerate,esp_lat,esp_jit,nesp,pesp,wn,dsp);
        sprate0 = sum(spt0(:,inx1:inx2),2);
        if rand > effect_probability   % in some cases there may no be an effect
            spt1 = spike_trains(mx,esp,baselinerate,esp_lat,esp_jit,nesp,pesp,wn,dsp);
        else
            spt1 = spike_trains(mx,esp,baselinerate,esp_lat,esp_jit,nesp,pesp*(1+effect_size),wn,dsp);
        end
        sprate1 = sum(spt1(:,inx1:inx2),2);
        [ROC(iC) P(iC) SE(iC)] = rocarea(sprate0,sprate1,'transform','scale','bootstrap',1000);
    end
    
    [p(iS) H(iS)] = signrank(ROC);
end

P = sum(H) / sno;   % statistical power

% -------------------------------------------------------------------------
function spt = spike_trains(mx,esp,baselinerate,esp_lat,esp_jit,nesp,pesp,wn,dsp)

% Simulate event train
E = esp:esp:mx;
evnt = length(E); % number of events

% Simulate spike train
S1 = randpoisson(mx/1000*baselinerate,mx);    % generate background Poisson spiking
S2 = repmat(E,nesp,1) + esp_lat + randn(nesp,length(E)) * esp_jit;  % generate 'evoked' spikes
S2 = sort(S2(:))';
prob_resp = rand(size(S2)) < pesp;
S2 = S2(prob_resp);   % probabilistic response
S = sort([S1 S2]);   % all spikes

% Curve fitting - parameters
a = 30;
c = baselinerate;
m = esp_lat;
s = 0.3;
s = fitoptions('Method','NonlinearLeastSquares',...   % non-linear least squares
    'Lower',[0,0,-Inf,-Inf],...    % lower bounds
    'Upper',[Inf,Inf,Inf,Inf],...   % upper bounds
    'Startpoint',[a,c,m,s],...   % initial values
    'Robust','on');    % robust fit
f = getmodel(s);

% Raster plot, PSTH
dt = 1;   % time resolution of bin raster (ms)
spt = raster_psth(E,S,wn,dt,evnt,esp,mx,f,s,dsp);

% -------------------------------------------------------------------------
function spt = raster_psth(Ev,S,wn,dt,evnt,esp,mx,f,s,dsp)

% Pre-align spikes
PerEvSpiks = struct('spikes',[],'event',[]);
SNo = evnt;    % SNo = evnt for full raster showing all events
for n = 1:SNo
    wnd = S(S>=(Ev(n)+wn(1))&S<=(Ev(n)+wn(2)));
    zeroed = wnd - Ev(n);   % reset values of all windows around 0 for plotting
    dimension = numel(zeroed);
    yupsilon = repmat(n,1,dimension);
    PerEvSpiks(n).spikes = zeroed;   % chunks of spikes + and - X ms around the event
    PerEvSpiks(n).event = yupsilon;
end

% Create scatter data of the spikes and PETH
ultravectorX = horzcat(PerEvSpiks(:).spikes);
ultravectorY = horzcat(PerEvSpiks(:).event);
time = wn(1):dt:wn(2);
spt = zeros(SNo,length(time));   % bin raster
ruX = round(ultravectorX-wn(1));
ruY = ultravectorY;
outinx = ruX==0;
ruX(outinx) = [];
ruY(outinx) = [];
spt(sub2ind([SNo length(time)],ruY,ruX)) = 1;
if dsp
    H1 = figure;   %#ok<NASGU> % raster plot and PSTH
    A1 = subplot(2,1,1);
    scatter(ultravectorX,ultravectorY,'.' , 'k')
    axis([wn(1) wn(2) 0 SNo]);
    hold on
    RR = plot([0 0],[0 SNo],'Color',[0.8 0 0]);
    RR(1).LineWidth = 2;
    hold off
    A2 = subplot(2,1,2);
    hist(ultravectorX,wn(2)-wn(1)+1) % alternatively, esp or t_window*2+1 for window size
    yl = ylim;
    hold on
    rl = plot([0 0], [0 yl(2)],'Color', [0.8 0 0]);
    rl(1).LineWidth = 2;
    hold off
    linkaxes([A1 A2],'x')
    
    % Plot events and spikes in time
    pse = zeros(1,mx+1000);   % pseudo-event train
    pse(ceil(Ev)) = 1;
    psu = zeros(1,mx+1000);   % pseudo-spike train
    psu(ceil(S)) = 1;
    H2 = figure; %#ok<NASGU>
    A1 = subplot(211);
    l = plot(pse,'Color',[0.8 0 0]);
    l(1).LineWidth = 2;
    ylim([-0.2 1.2])
    A2 = subplot(212);
    ll = plot(psu,'k');
    ll(1).LineWidth = 2;
    ylim([-0.2 1.2])
    linkaxes([A1 A2])
    box off
    set(gca,'xticklabel',[]);
    
    % Cross-correlation
    t_window = max(abs(wn));
    [R, lags] = xcorr(psu,pse,t_window);   % E -> U
    [fun , ~] = fit(lags',R',f,s);
    figure
    c = plot(lags,R,'Color',[0 0 0.8]);
    xlim(wn);
    c(1).LineWidth = 2;
    hold on
    L = plot(fun,'k');   % overlay the fitted curve
    L(1).LineWidth = 2;
    hold on
    yl = ylim;
    rl = plot([0 0],[0 yl(2)],'Color',[0.8 0 0]);
    rl(1).LineWidth = 2;
    hold off
    box off
    set(gca,'TickDir','out');
    ylabel('cross-correlation');
end

% -------------------------------------------------------------------------
function f = getmodel(s)

f = fittype('c+a*exp(-((x-m)/s)^2)','options',s);