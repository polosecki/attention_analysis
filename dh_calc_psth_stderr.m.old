function [psth] = dh_calc_psth_stderr(fid,bid,ts,te,clussel,tres,Num,Den,fdelay,dec)
% Calculate the PSTH for an opened DH5-file. Arguments are
% fid - handle to an open file
% bid - spike block id to calculate psth from
% ts, te - time intervals
% clussel - include only spikes with specified cluster number(s)
%           -1 means no cluster selection (all spikes)
% tres - time resolution (in nanoseconds)
% Num,Den - coefficients of FIR or IIR filter used
%           for histogram-like filtering
% fdelay - FIR filter delay to compensate for
%          (use zero for non linear-phase filters)
% dec - downsample the filtered and non-filtered result
%       by this factor

nspk = double(dhfun('GETNUMBERSPIKES',fid,bid));

tstamps = dhfun('READSPIKEINDEX',fid,bid,1,nspk);
% If there was provided a cluster selection, perform it
if clussel >= 0
    clusinfo = dhfun('READSPIKECLUSTER',fid,bid,1,nspk);

    % Select only spikes with given cluster numbers
    tstamps = tstamps(ismember(clusinfo,clussel));
end

% Find the longest trial
maxtl = max(te-ts);
% Find the number of bins in this trial
maxbw = ceil(maxtl/tres)+1;

% Prepare the resulting psth array and supplemental data
contrib = zeros(1,maxbw);
psthdata = zeros(1,maxbw);
psthsum_f = zeros(1,maxbw);
psthsumsq_f = zeros(1,maxbw);
%fprintf('maxbw = %d\n',maxbw);

% Process, trial by trial
for tk = 1:length(ts)
    tlen = te(tk)-ts(tk);
    tlenbw = ceil(tlen/tres)+1;
    odata = zeros(1,tlenbw);
    otstamps = tstamps(tstamps>=ts(tk) & tstamps <= te(tk));
    otstamps = otstamps - ts(tk); %Refer all spikes to the start of trial
    odata(floor(otstamps/tres) + 1) = 1;
    contrib(1:tlenbw) = contrib(1:tlenbw) + 1;
    %fprintf('tk= %d tlenbw = %d lodata=%d\n',tk,tlenbw,length(odata));
    psthdata(1:tlenbw) = psthdata(1:tlenbw) + odata;

    % Calculate the filtered PSTH
    [odata_f,zf] = filter(Num,Den,odata);
    % Compensate for the filter delay
    odata_f = [odata_f(fdelay+1:end) filter(Num,Den,zeros(1,fdelay),zf)];
    psthsum_f(1:tlenbw) = psthsum_f(1:tlenbw) + odata_f;
    psthsumsq_f(1:tlenbw) = psthsumsq_f(1:tlenbw) + odata_f.^2;
end

% Calculate the standard deviation
coco = contrib>1;
icoco = contrib==1;
dvar = zeros(size(psthsum_f));
dvar(coco) = (psthsumsq_f(coco) - (psthsum_f(coco).^2)./contrib(coco))./(contrib(coco)-1);
psth_stddev = sqrt(dvar);
psth_stderr = zeros(size(psth_stddev));
psth_stderr(coco) = psth_stddev(coco) ./ sqrt(contrib(coco)-1);
psth_stderr(icoco) = max(psth_stderr);
% Apply the contribution info
psthdata = psthdata ./ contrib;

% Divide the values by time resolution, thereby
% obtaining the actual density of spikes per second
% (Firing rate)
psthdata = psthdata ./ tres .* 1e9; % tres is specified in nanoseconds
% Do the same with stddev and stderr
psth_stddev = psth_stddev ./ tres .* 1e9;
psth_stderr = psth_stderr ./ tres .* 1e9;

% Calculate the filtered PSTH
[psthdata_f,zf] = filter(Num,Den,psthdata);
% Compensate for the filter delay
adata = filter(Num,Den,zeros(1,fdelay),zf);
psthdata_f = [psthdata_f(fdelay+1:end) adata];
% Decimate the filtered psth data and contribution
psthdata_f = psthdata_f(1:dec:end);
psth_stddev = psth_stddev(1:dec:end);
psth_stderr = psth_stderr(1:dec:end);
contrib = contrib(1:dec:end);
% Calculate the non-filtered, decimated PSTH
if dec > 1
    Numd = ones(1,dec)/dec;
    fdelayd = floor(dec/2);
    [psthdata_nf,zf] = filter(Numd,1,psthdata);
    adata = filter(Numd,1,zeros(1,fdelayd),zf);
    psthdata_nf = [psthdata_nf(fdelayd+1:end) adata];
    psthdata_nf = psthdata_nf(1:dec:end);
else
    psthdata_nf = psthdata;
end

psth.data = psthdata_f;
psth.data_nf = psthdata_nf;
psth.contrib = contrib;
psth.stddev = psth_stddev;
psth.stderr = psth_stderr;
