%% 0. Setup (repo-relative paths)
clear; clc; close all;
repoDir  = fileparts(mfilename('fullpath'));
miDir    = fullfile(repoDir,'MI_data_scripts','MI_data_scripts');
biosigDir= fullfile(repoDir,'biosig');
dataDir  = fullfile(miDir,'Subject_006_Session_006_TESS_Online_Visual');
assert(isfolder(biosigDir),'biosig not found: %s',biosigDir);
assert(isfolder(dataDir),  'Data folder not found: %s',dataDir);
addpath(genpath(biosigDir));
addpath(miDir);                      % sload.m, topoplot_.m, selectedChannels.mat
files = dir(fullfile(dataDir,'*.gdf'));
[~,ord] = sort({files.name}); files = files(ord);
nRuns = numel(files);
assert(nRuns==4,'Expected 4 runs, found %d',nRuns);
which sload topoplot_                % confirm the intended copies are used

%% 1. Parameters
order = 4; band = [8 12]; keepCh = 1:32;
code.cue   = struct('RH',769,  'LH',770);
code.start = struct('RH',7691, 'LH',7701);
code.miss  = struct('RH',7692, 'LH',7702);
code.hit   = struct('RH',7693, 'LH',7703);
lastWin_s = 0.5;


%% 2. Channel labels, Laplacian neighbors, spatial filter matrices
chanLabels = {'FP1','FPZ','FP2','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6', ...
              'M1','T7','C3','CZ','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3', ...
              'PZ','P4','P8','POZ','O1','OZ','O2'};
assert(numel(chanLabels) == 32);

nbrNames = struct( ...
 'FP1', {{'FPZ','F3','F7'}}, ...
 'FPZ', {{'FP1','FP2','FZ'}}, ...
 'FP2', {{'FPZ','F4','F8'}}, ...
 'F7',  {{'FP1','F3','FC5'}}, ...
 'F3',  {{'FP1','F7','FZ','FC5','FC1'}}, ...
 'FZ',  {{'FPZ','F3','F4','FC1','FC2'}}, ...
 'F4',  {{'FP2','FZ','F8','FC2','FC6'}}, ...
 'F8',  {{'FP2','F4','FC6'}}, ...
 'FC5', {{'F7','F3','FC1','T7','C3'}}, ...
 'FC1', {{'F3','FZ','FC5','FC2','C3','CZ'}}, ...
 'FC2', {{'FZ','F4','FC1','FC6','CZ','C4'}}, ...
 'FC6', {{'F4','F8','FC2','C4','T8'}}, ...
 'M1',  {{'T7','CP5'}}, ...
 'T7',  {{'FC5','C3','CP5'}}, ...
 'C3',  {{'FC5','FC1','T7','CZ','CP5','CP1'}}, ...
 'CZ',  {{'FC1','FC2','C3','C4','CP1','CP2'}}, ...
 'C4',  {{'FC2','FC6','CZ','T8','CP2','CP6'}}, ...
 'T8',  {{'FC6','C4','CP6'}}, ...
 'M2',  {{'T8','CP6'}}, ...
 'CP5', {{'T7','C3','CP1','P7','P3'}}, ...
 'CP1', {{'C3','CZ','CP5','CP2','P3','PZ'}}, ...
 'CP2', {{'CZ','C4','CP1','CP6','PZ','P4'}}, ...
 'CP6', {{'C4','T8','CP2','P4','P8'}}, ...
 'P7',  {{'CP5','P3','O1'}}, ...
 'P3',  {{'CP5','CP1','P7','PZ','O1','POZ'}}, ...
 'PZ',  {{'CP1','CP2','P3','P4','POZ'}}, ...
 'P4',  {{'CP2','CP6','PZ','P8','POZ','O2'}}, ...
 'P8',  {{'CP6','P4','O2'}}, ...
 'POZ', {{'PZ','P3','P4','OZ'}}, ...
 'O1',  {{'P7','P3','POZ','OZ'}}, ...
 'OZ',  {{'POZ','O1','O2'}}, ...
 'O2',  {{'P8','P4','POZ','OZ'}});

nbr = cell(32,1);
for i = 1:32
    names = nbrNames.(upper(chanLabels{i}));
    nbr{i} = cellfun(@(n) find(strcmpi(chanLabels, n)), names);
end
assert(all(cellfun(@numel, nbr) >= 1), 'A channel has no neighbors');
assert(~any(arrayfun(@(i) ismember(i, nbr{i}), 1:32)), 'A channel lists itself');

W_car = eye(32) - ones(32)/32;
W_lap = eye(32);
for i = 1:32
    W_lap(i, nbr{i}) = -1/numel(nbr{i});
end
fprintf('\n=== 2. Spatial filters defined ===\n');
fprintf('Neighbors per channel:\n');
for i = 1:32
    fprintf('  %-4s <- %s\n', chanLabels{i}, strjoin(chanLabels(nbr{i}), ', '));
end
%% 3. Load and filter every run
EEG = cell(nRuns,1);  HDR = cell(nRuns,1);
MU  = cell(nRuns,1);  CAR = cell(nRuns,1);  LAP = cell(nRuns,1);
for r = 1:nRuns
    fn = fullfile(dataDir, files(r).name);
    [s, h] = sload(fn);
    fs_r = h.SampleRate;
    if r == 1, fs = fs_r; chanLabels = h.Label(keepCh); end

    assert(fs_r == fs, 'Sampling rate differs in run %d', r);
    assert(isequal(h.Label(keepCh), chanLabels), 'Channel labels differ in run %d', r);
    assert(size(s,2) >= 32, 'Run %d has fewer than 32 columns', r);

    eeg = s(:, keepCh);
    assert(~any(isnan(eeg(:))), 'NaN in raw data, run %d', r);

    [b, a] = butter(order, band/(fs/2), 'bandpass');
    mu = filtfilt(b, a, eeg);

    EEG{r} = eeg;  HDR{r} = h;
    MU{r}  = mu;
    CAR{r} = mu * W_car';
    LAP{r} = mu * W_lap';

    fprintf('\n--- Run %d: %s ---\n', r, files(r).name);
    fprintf('Raw size %d x %d -> kept %d x %d | fs %g Hz | length %.1f s\n', ...
        size(s,1), size(s,2), size(eeg,1), size(eeg,2), fs_r, size(eeg,1)/fs_r);
    fprintf('Dropped columns: %s\n', strjoin(h.Label(33:end), ', '));
    fprintf('Filter stable: %d | NaN/Inf after filter: %d\n', ...
        all(abs(roots(a)) < 1), any(~isfinite(mu(:))));
    fprintf('Ch1 (%s) mu-band power fraction: before %.3f, after %.3f\n', chanLabels{1}, ...
        bandpower(eeg(:,1), fs, band)/bandpower(eeg(:,1), fs, [0.5 fs/2-1]), ...
        bandpower(mu(:,1),  fs, band)/bandpower(mu(:,1),  fs, [0.5 fs/2-1]));
    fprintf('CAR max |mean across ch|: %.2e | Laplacian NaN: %d\n', ...
        max(abs(mean(CAR{r},2))), any(isnan(LAP{r}(:))));
end

% Check that the matrix form of the Laplacian equals the explicit loop (run 1)
lap_loop = MU{1};
for i = 1:32, lap_loop(:,i) = MU{1}(:,i) - mean(MU{1}(:, nbr{i}), 2); end
fprintf('\nMatrix vs loop Laplacian max diff (run 1): %.2e\n', max(abs(lap_loop(:) - LAP{1}(:))));

%% 4. Inspect events in every run
fprintf('\n=== 4. Events per run ===\n');
for r = 1:nRuns
    typ = HDR{r}.EVENT.TYP;  pos = HDR{r}.EVENT.POS;
    [u,~,k] = unique(typ);  cnt = accumarray(k,1);
    fprintf('Run %d: %d events\n', r, numel(typ));
    for i = 1:numel(u)
        fprintf('   code 0x%04X (%d): %d events\n', u(i), u(i), cnt(i));
    end
    n = min(12, numel(typ));
    fprintf('   First %d events (code, time s):\n', n);
    disp([typ(1:n), pos(1:n)/fs]);
    if isfield(HDR{r}.EVENT,'DUR') && ~isempty(HDR{r}.EVENT.DUR)
        fprintf('   First durations (samples): %s\n', mat2str(HDR{r}.EVENT.DUR(1:n)'));
    end
end

%% 5. Count trials (cue-based) per class and run
cls = {'RH','LH'};
nTr = zeros(nRuns,2);
for r = 1:nRuns
    typ = HDR{r}.EVENT.TYP;
    for c = 1:2, nTr(r,c) = sum(typ==code.cue.(cls{c})); end
end
disp(array2table(nTr,'VariableNames',cls,'RowNames',compose('run%d',1:nRuns)));
fprintf('Total RH %d, LH %d, all %d\n',sum(nTr(:,1)),sum(nTr(:,2)),sum(nTr(:)));



%% 6. Extract task trials (feedback start -> hit/miss), variable length
T = struct();
for c = 1:2
    k = cls{c};
    T.(k).mu={}; T.(k).car={}; T.(k).lap={}; T.(k).run=[]; T.(k).outcome=[]; T.(k).len_s=[];
    for r = 1:nRuns
        typ=HDR{r}.EVENT.TYP(:); pos=HDR{r}.EVENT.POS(:);
        iS = find(typ==code.start.(k));
        for t = 1:numel(iS)
            j = iS(t)+1;                         % next event = end marker
            if j>numel(typ) || ~ismember(typ(j),[code.miss.(k) code.hit.(k)]), continue; end
            idx = pos(iS(t)):pos(j);
            T.(k).mu{end+1}=MU{r}(idx,:);  T.(k).car{end+1}=CAR{r}(idx,:);
            T.(k).lap{end+1}=LAP{r}(idx,:);
            T.(k).run(end+1)=r;  T.(k).outcome(end+1)=(typ(j)==code.hit.(k));
            T.(k).len_s(end+1)=numel(idx)/fs;
        end
    end
    fprintf('%s: %d trials, hit %d, miss %d, length %.2f-%.2f s\n',k,numel(T.(k).run), ...
        sum(T.(k).outcome),sum(~T.(k).outcome),min(T.(k).len_s),max(T.(k).len_s));
end
% Q: trials per class extracted must equal cue counts in section 5, else investigate


%% 7. Per-trial mu power (mean of squared samples), raw vs CAR
cls = {'RH','LH'};  filt = {'mu','car'};  filtName = {'No spatial filter','CAR'};
nL = round(lastWin_s*fs);
P = struct();  GA = struct();
for c = 1:2
    k = cls{c};  nT = numel(T.(k).run);
    for f = 1:2
        X = T.(k).(filt{f});
        Pfull = zeros(32,nT);  Plast = nan(32,nT);
        for t = 1:nT
            Pfull(:,t) = mean(X{t}.^2, 1)';                  % whole task, 32 x 1
            if size(X{t},1) >= nL
                Plast(:,t) = mean(X{t}(end-nL+1:end,:).^2, 1)';   % last 0.5 s
            end
        end
        P.(k).(filt{f})  = Pfull;                            % [channels x trials]
        GA.(k).(filt{f}) = mean(Plast, 2, 'omitnan');        % [channels x 1]
        assert(~any(isnan(Plast(:))), 'A %s trial was shorter than %.1f s', k, lastWin_s);
        fprintf('%s %-3s: per-trial power %s | GA range %.3g to %.3g\n', k, filt{f}, ...
            mat2str(size(Pfull)), min(GA.(k).(filt{f})), max(GA.(k).(filt{f})));
    end
end
% Units: signal units squared. Check h.PhysDim (likely uV, so uV^2)
disp(HDR{1}.PhysDim);     % amplitude units for the power values (likely uV, so uV^2)

%% 8. Topoplots of last-0.5 s grand-average mu power
assert(exist('GA','var')==1 && isfield(GA,'RH'), 'Run sections 0-7 first (GA missing)');
S = load(fullfile(miDir,'selectedChannels.mat'));
selCh = S.selectedChannels;
assert(numel(selCh)==32, 'Expected 32 channel locations, got %d', numel(selCh));
assert(isequal(upper({selCh.labels}), upper(chanLabels(:)')), 'Channel order mismatch');
fprintf('selCh urchan (original montage indices, informational): %s\n', mat2str([selCh.urchan]));

cls = {'RH','LH'};  fl = {'mu','car'};  flName = {'No spatial filter','CAR'};
figure('Name','Mu power, last 0.5 s','Position',[100 100 900 700]);
for f = 1:2
    lo = min([GA.RH.(fl{f}); GA.LH.(fl{f})]);
    hi = max([GA.RH.(fl{f}); GA.LH.(fl{f})]);
    for c = 1:2
        subplot(2,2,(c-1)*2+f);
        topoplot_(GA.(cls{c}).(fl{f}), selCh, 'maplimits', [lo hi]);
        colorbar;
        title(sprintf('%s | %s', cls{c}, flName{f}));
    end
end
sgtitle('Grand-average mu power (8-12 Hz), last 0.5 s of task, n=40 per class');