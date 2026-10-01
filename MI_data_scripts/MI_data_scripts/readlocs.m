function [eloc, labels, Th, Rd, indices] = readlocs(chans, varargin)
% Minimal stand-in for EEGLAB readlocs, struct input only.
    assert(isstruct(chans) && all(isfield(chans,{'labels','theta','radius'})), ...
        'readlocs stub needs a struct with labels, theta, radius');
    eloc    = chans;
    labels  = {chans.labels};
    hasLoc  = ~cellfun(@isempty,{chans.theta}) & ~cellfun(@isempty,{chans.radius});
    indices = find(hasLoc);
    Th      = double([chans(indices).theta]);     % degrees, topoplot_ converts to radians
    Rd      = double([chans(indices).radius]);
end