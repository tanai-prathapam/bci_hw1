%% 0. Setup (repo-relative paths, same convention as hw1_q1.m)
clear; clc; close all;
rng(1);
repoDir = fileparts(mfilename('fullpath'));
errDir  = fullfile(repoDir,'ErrP_data_scripts','ErrP_data_scripts');
miDir   = fullfile(repoDir,'MI_data_scripts','MI_data_scripts');   % fallback topoplot_ only
outDir  = fullfile(repoDir,'output_q2');
assert(isfolder(errDir),'ErrP folder not found: %s',errDir);
if ~isfolder(outDir), mkdir(outDir); end
addpath(errDir);
which topoplot
repoCommitRead = 'f0c2085e59ec6c4d9dbb953bc0101a835ad43e10';  % main HEAD when inputs were read
fprintf('Repo commit read: %s\nMATLAB %s\n', repoCommitRead, version);

%% 1. Parameters
winSec     = [-0.2 0.8];     % analysis window (s) per assignment
baseSec    = [-0.2 0];       % baseline window (s), [start, end)
doBaseline = true;
ernWin     = [0.05 0.40];    % PROVISIONAL ERN search window (s): confirm on plotted waveform
peEnd      = 0.80;           % Pe searched from ERN to this time (s)
ernRule    = 'min';          % 'min' = most negative prominent minimum, 'first' = first prominent minimum
promFrac   = 0.10;           % min prominence as a fraction of post-0 peak-to-peak range
ccaOnCAR   = false;          % CCA on baseline-corrected non-CAR data by default
nComp      = 5;
tol        = 1e-9;

%% 2. Load and validate
D = load(fullfile(errDir,'ErrP_data_HW1.mat'));
assert(isfield(D,'trainingEpochs'),'Expected trainingEpochs');
E = D.trainingEpochs;
X = E.rotation_data;  label = E.label(:);  mag = E.magnitude(:);
C = load(fullfile(errDir,'ErrP_channels.mat'));
P = C.params;  CL = P.chanlocs;  fs = P.fsamp;  t = P.epochTime(:);
[nT,nCh,nTr] = size(X);
assert(isequal([nT nCh nTr],[1024 32 320]),'Unexpected size %s',mat2str([nT nCh nTr]));
assert(nCh == numel(CL),'Channels (%d) ~= chanlocs (%d)',nCh,numel(CL));
assert(isequal(P.eegChannels(:)',1:nCh),'Not all 32 channels are EEG');
assert(numel(t)==nT && abs(median(diff(t))-1/fs)<tol,'Time axis mismatch');
assert(abs(t(P.epochOnset))<tol,'epochOnset is not t = 0');
assert(~any(isnan(X(:))),'NaN in data');
assert(numel(label)==nTr && numel(mag)==nTr,'label/magnitude length mismatch');
assert(all(ismember(label,[0 1])),'Unexpected labels');
assert(all(mag(label==0)==0) && all(mag(label==1)>0),'label/magnitude inconsistent');
chanNames = strtrim(arrayfun(@(c) char(c.labels),CL,'UniformOutput',false));
czIdx = find(strcmpi(chanNames,'Cz'));  fczIdx = find(strcmpi(chanNames,'FCz'));
assert(isscalar(czIdx) && isscalar(fczIdx),'Cz/FCz lookup failed');
fprintf('fs %g Hz | Cz idx %d | FCz idx %d | t(%d)=%.4f s | t(1)=%.4f, t(end)=%.4f\n', ...
    fs,czIdx,fczIdx,P.epochOnset,t(P.epochOnset),t(1),t(end));
fprintf('Units: not stated in file; |x| median/99%%/99.9%% = 2.62/12.38/17.83 -> plotting as a.u.\n');
disp('--- params.spectralFilter / spatialFilter / resample (check for prior CAR or other filtering) ---');
disp(P.spectralFilter); disp(P.spatialFilter); disp(P.resample);

%% 3. Window and baseline
win     = t >= winSec(1)-tol & t <= winSec(2)+tol;
tw      = t(win);  nWin = nnz(win);  winIdx = find(win);
baseIdx = t >= baseSec(1)-tol & t < baseSec(2)-tol;
fprintf('Window [%.3f, %.3f] s: %d samples (first %.4f s, last %.4f s)\n',winSec(1),winSec(2),nWin,tw(1),tw(end));
fprintf('Baseline [%.2f, %.2f) s: %d samples | doBaseline = %d\n',baseSec(1),baseSec(2),nnz(baseIdx),doBaseline);
Xb = X;
if doBaseline, Xb = X - mean(X(baseIdx,:,:),1); end

%% 4. CAR (computed on epochs before averaging)
Xc = Xb - mean(Xb,2);
carMaxMean = max(abs(mean(Xc(win,:,:),2)),[],'all');
assert(carMaxMean < 1e-9,'CAR channel mean not ~0');
fprintf('CAR: max |mean across %d channels| = %.2e\n',nCh,carMaxMean);
dsName  = {'raw','car'};  dsTitle = {'No spatial filter','CAR'};  dsData = {Xb,Xc};

%% 5. Per-magnitude Cz grand average, ERN/Pe detection, figures
magList = unique(mag)';  errIdx = label==1;
groups = {'All error', errIdx};
for m = magList(magList>0)
    groups(end+1,:) = {sprintf('%d deg error',m), mag==m}; %#ok<SAGROW>
end
rows = {};  DET = struct();
for f = 1:2
    Xf = dsData{f};
    fig = figure('Name',['Cz GA | ' dsTitle{f}],'Position',[100 100 950 600]); hold on;
    hs = []; leg = {};
    for m = magList
        idx = mag==m;
        ga  = mean(squeeze(Xf(win,czIdx,idx)),2);
        hs(end+1) = plot(tw*1000,ga,'LineWidth',1.2); %#ok<SAGROW>
        if m==0, leg{end+1} = sprintf('0 deg (correct) (n=%d)',nnz(idx)); %#ok<SAGROW>
        else,    leg{end+1} = sprintf('%d deg (n=%d)',m,nnz(idx)); end %#ok<SAGROW>
    end
    gaErr = mean(squeeze(Xf(win,czIdx,errIdx)),2);
    hs(end+1) = plot(tw*1000,gaErr,'k--','LineWidth',2);
    leg{end+1} = sprintf('All error (n=%d)',nnz(errIdx));
    for g = 1:size(groups,1)
        w = mean(squeeze(Xf(win,czIdx,groups{g,2})),2);
        d = detectErrP(w,tw,ernWin,ernRule,promFrac,peEnd);
        if g==1, DET.(dsName{f}) = d; end
        rows(end+1,:) = {dsTitle{f},groups{g,1},nnz(groups{g,2}),d.tERN*1000,d.aERN,d.tPe*1000,d.aPe}; %#ok<SAGROW>
    end
    d = DET.(dsName{f});
    yl = ylim;
    patch([ernWin fliplr(ernWin)]*1000,[yl(1) yl(1) yl(2) yl(2)],[.85 .85 .85], ...
        'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
    uistack(findobj(gca,'Type','patch'),'bottom');
    if ~isnan(d.tERN), plot(d.tERN*1000,d.aERN,'kv','MarkerFaceColor','r','MarkerSize',9,'HandleVisibility','off'); text(d.tERN*1000,d.aERN,sprintf('  ERN %.0f ms',d.tERN*1000)); end
    if ~isnan(d.tPe),  plot(d.tPe*1000,d.aPe,'k^','MarkerFaceColor','g','MarkerSize',9,'HandleVisibility','off');   text(d.tPe*1000,d.aPe,sprintf('  Pe %.0f ms',d.tPe*1000)); end
    xline(0,'k:','HandleVisibility','off');  yline(0,'k:','HandleVisibility','off');
    xlabel('Time from trigger (ms)'); ylabel('Amplitude at Cz (a.u., unit unconfirmed)');
    title(sprintf('Cz grand average by rotation magnitude | %s | baseline %d',dsTitle{f},doBaseline));
    legend(hs,leg,'Location','best'); grid on;
    saveFig(fig,outDir,sprintf('q2_cz_ga_%s',dsName{f}));
end
LAT = cell2table(rows,'VariableNames',{'Filter','Group','nTrials','ERN_ms','ERN_amp','Pe_ms','Pe_amp'});
disp(LAT);
fprintf('ERN shift after CAR: %+.2f ms | Pe shift after CAR: %+.2f ms (all-error GA, Cz)\n', ...
    (DET.car.tERN-DET.raw.tERN)*1000,(DET.car.tPe-DET.raw.tPe)*1000);

%% 6. Topoplots at ERN and Pe (all error trials), raw and CAR
for f = 1:2
    d = DET.(dsName{f});
    assert(~isnan(d.tERN) && ~isnan(d.tPe),'ERN/Pe not found for %s; adjust ernWin/promFrac',dsName{f});
    iE = winIdx(d.iERN);  iP = winIdx(d.iPe);
    vE = squeeze(mean(dsData{f}(iE,:,errIdx),3)); vE = vE(:);
    vP = squeeze(mean(dsData{f}(iP,:,errIdx),3)); vP = vP(:);
    assert(numel(vE)==nCh && numel(vP)==nCh);
    L = max(abs([vE;vP]));                       % symmetric shared limits [-L L]
    fprintf('%s: ERN %.1f ms (spatial mean %.3g), Pe %.1f ms (spatial mean %.3g), color limits +/-%.3g\n', ...
        dsTitle{f},d.tERN*1000,mean(vE),d.tPe*1000,mean(vP),L);
    fig = figure('Name',['Topo | ' dsTitle{f}],'Position',[100 100 900 420]);
    subplot(1,2,1); safeTopo(vE,CL,[-L L],miDir); cb = colorbar; ylabel(cb,'Amplitude (a.u.)');
    title(sprintf('ERN %.0f ms | %s',d.tERN*1000,dsTitle{f}));
    subplot(1,2,2); safeTopo(vP,CL,[-L L],miDir); cb = colorbar; ylabel(cb,'Amplitude (a.u.)');
    title(sprintf('Pe %.0f ms | %s',d.tPe*1000,dsTitle{f}));
    sgtitle(sprintf('Grand average, all error trials (n=%d), symmetric shared scale',nnz(errIdx)));
    saveFig(fig,outDir,sprintf('q2_topo_ern_pe_%s',dsName{f}));
end

%% 7. CCA spatial filter (error vs correct templates)
% NOTE (leakage): descriptive only. Templates and filters are fit on the same
% trials that are plotted, so canonical correlations are optimistic. No held-out test.
Xs = Xb; if ccaOnCAR, Xs = Xc; end
tplErr = mean(Xs(win,:,label==1),3);   % nWin x nCh
tplCor = mean(Xs(win,:,label==0),3);
concat_data = zeros(nWin*nTr,nCh);  concat_ga = zeros(nWin*nTr,nCh);
for k = 1:nTr
    r = (k-1)*nWin + (1:nWin);
    concat_data(r,:) = Xs(win,:,k);
    if label(k)==1, concat_ga(r,:) = tplErr; else, concat_ga(r,:) = tplCor; end
end
assert(~any(isnan(concat_data(:))) && ~any(isnan(concat_ga(:))));
fprintf('CCA: %d trials x %d samples, rank(data)=%d, rank(template)=%d, ccaOnCAR=%d\n', ...
    nTr,nWin,rank(concat_data),rank(concat_ga),ccaOnCAR);
[A,~,rho] = canoncorr(concat_data,concat_ga);     % A: channels x components
assert(size(A,1)==nCh && size(A,2)>=nComp,'Fewer than %d CCA components (rank deficient)',nComp);
A = A(:,1:nComp);  rho = rho(1:nComp);
for k = 1:nComp, if A(fczIdx,k) < 0, A(:,k) = -A(:,k); end, end   % sign: FCz weight positive
fprintf('Canonical correlations 1-%d: %s\n',nComp,mat2str(rho,4));
fig = figure('Name','CCA weights','Position',[100 100 1100 650]);
for k = 1:nComp
    subplot(2,3,k); Lk = max(abs(A(:,k)));
    safeTopo(A(:,k),CL,[-Lk Lk],miDir); cb = colorbar; ylabel(cb,'Filter weight (a.u.)');
    title(sprintf('CCA %d (r = %.3f)',k,rho(k)));
end
sgtitle('CCA spatial filter weights (error vs correct templates), per-component symmetric scale, FCz sign-normalized');
saveFig(fig,outDir,'q2_cca_weights_all');
for k = 1:nComp
    fig = figure('Position',[100 100 480 420]); Lk = max(abs(A(:,k)));
    safeTopo(A(:,k),CL,[-Lk Lk],miDir); cb = colorbar; ylabel(cb,'Filter weight (a.u.)');
    title(sprintf('CCA %d (r = %.3f)',k,rho(k)));
    saveFig(fig,outDir,sprintf('q2_cca_comp%d',k)); close(fig);
end

%% 8. Save results
cnt = array2table([sum(label==0 & mag==0), arrayfun(@(m) sum(label==1 & mag==m),magList(magList>0))], ...
    'VariableNames',["correct_0deg" compose("error_%ddeg",magList(magList>0))]);
CCA = table((1:nComp)',rho(:),'VariableNames',{'Component','CanonicalCorr'});
R = struct('commit',repoCommitRead,'matlab',version,'fs',fs,'window_s',winSec,'nWinSamples',nWin, ...
    'doBaseline',doBaseline,'baseSec',baseSec,'ernWin',ernWin,'ernRule',ernRule,'promFrac',promFrac, ...
    'latencies',LAT,'counts',cnt,'cca_rho',rho,'cca_weights',A,'chanNames',{chanNames});
save(fullfile(outDir,'results.mat'),'R');
writetable(LAT,fullfile(outDir,'latencies.csv'));
writetable(cnt,fullfile(outDir,'trial_counts.csv'));
writetable(CCA,fullfile(outDir,'cca_corr.csv'));
disp(cnt); disp(CCA);
fprintf('Saved outputs to %s\n',outDir);

%% Local functions
function d = detectErrP(w,tw,ernWin,ernRule,promFrac,peEnd)
    w = w(:);
    d = struct('iERN',NaN,'tERN',NaN,'aERN',NaN,'iPe',NaN,'tPe',NaN,'aPe',NaN);
    post = tw > 0;
    prom = promFrac*(max(w(post))-min(w(post)));
    iw   = find(tw >= ernWin(1) & tw <= ernWin(2));
    [pk,loc] = findpeaks(-w(iw),'MinPeakProminence',prom);
    if isempty(loc), return; end
    if strcmp(ernRule,'first'), k = 1; else, [~,k] = max(pk); end
    iN = iw(loc(k));
    ia = find(tw > tw(iN) & tw <= peEnd);
    if numel(ia) >= 3
        [~,loc2] = findpeaks(w(ia),'MinPeakProminence',prom);
    else, loc2 = []; end
    d.iERN = iN; d.tERN = tw(iN); d.aERN = w(iN);
    if ~isempty(loc2)
        iP = ia(loc2(1)); d.iPe = iP; d.tPe = tw(iP); d.aPe = w(iP);
    end
end

function safeTopo(v,CL,lim,miDir)
    try
        topoplot(v,CL,'maplimits',lim);
    catch ME
        warning('topoplot failed (%s); falling back to topoplot_',ME.message);
        addpath(miDir); topoplot_(v,CL,'maplimits',lim);
    end
end

function saveFig(fig,outDir,name)
    exportgraphics(fig,fullfile(outDir,[name '.png']),'Resolution',300);
end