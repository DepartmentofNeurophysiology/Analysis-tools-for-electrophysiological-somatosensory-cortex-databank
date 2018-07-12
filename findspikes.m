function [spiketimes, spikeindices] = findspikes(dt, v, threshold, timestart)
% Find spikes from membrane potential trace

% INPUT
% dt: binsize of recordings, used for the calculation of spike times (so
% determines unit of spike times)
% v: vector with membrane potential recording
% threshold: membrane potential threshold above which to consider a spike
% timestart (optional): start time recording (if not zero): should have
% same unit as dt

% OUTPUT
% spiketimes: Nspikes x 

if nargin == 3
    timestart = 0;
end




    % parts of v above threshold
    thr=find(v>threshold);
    if isempty(thr)
        spiketimes = [];
        spikeindices = [];
        return
    end

    % find jumps: which belongs to one spike
    for i=2:length(thr)
        b(i-1)=thr(i)-thr(i-1);
    end
    jumps=find(b>1);

    if ~isempty(jumps)
        % so to one spike belong thrsom(jumpss(i-1)+1) 'till thrsom(jumpss(i)). Where
        % is the spiketime? Middle or max...

        spikeindices=zeros(1,length(jumps)+1);
        maxv=zeros(1,length(jumps)+1); 

        %middle
        %spiketimes(1)=dt*((thr(jumps(1)) - thr(1))/2 + thr(1)-1);
        %for i=2:length(jumps);
        %    spiketimes(i)=dt*((thr(jumps(i-1)+1)+(thr(jumps(i))-thr(jumps(i-1)+1))/2)-1);
        %end
        %spiketimes(length(jumps)+1)=dt*((thr(jumps(length(jumps))+1)+(thr(end)-thr(jumps(length(jumps))+1))/2)-1);

        % max
        maxv(1)=max(v(thr(1):thr(jumps(1))));
        if length(find(v(thr(1):thr(jumps(1)))==maxv(1))+thr(1)-1) == 1
            spikeindices(1)=find(v(thr(1):thr(jumps(1)))==maxv(1))+thr(1)-1;
        else
            indtemp = find(v(thr(1):thr(jumps(1)))==maxv(1))+thr(1)-1;
            indtemp = indtemp(1);
            spikeindices(1) = indtemp;
        end

        for i=2:length(jumps)
            maxv(i)=max(v(thr(jumps(i-1)+1):thr(jumps(i))));
            if length(find(v(thr(jumps(i-1)+1):thr(jumps(i)))== maxv(i))+thr(jumps(i-1)+1)-1)==1
                spikeindices(i)=find(v(thr(jumps(i-1)+1):thr(jumps(i)))== maxv(i))+thr(jumps(i-1)+1)-1;
            else %more than one value: take the first
                indtemp = find(v(thr(jumps(i-1)+1):thr(jumps(i)))== maxv(i))+thr(jumps(i-1)+1)-1;
                indtemp = indtemp(1);
                spikeindices(i)= indtemp;
            end
        end
        maxv(length(jumps)+1)=max(v(thr(jumps(end)+1):thr(end)));
        if length(find(v(thr(jumps(end)+1):thr(end))==maxv(length(jumps)+1))+thr(jumps(end)+1)-1) == 1
            spikeindices(length(jumps)+1)=find(v(thr(jumps(end)+1):thr(end))==maxv(length(jumps)+1))+thr(jumps(end)+1)-1;
        else
            indtemp = find(v(thr(jumps(end)+1):thr(end))==maxv(length(jumps)+1))+thr(jumps(end)+1)-1;
            indtemp = indtemp(1);
            spikeindices(length(jumps)+1) = indtemp;
        end
        
        spiketimes=(spikeindices-1)*dt+timestart;
    else
        spiketimes = [];
        spikeindices = [];
    end
end