function [oax,hg] = psth_plot(psth,tres,dopts,gcolor,contrib_color,ax,elec,clust,ylab)

tres = tres / 1e9; % Make tres be in seconds
ld = length(psth.data);
t = 0:tres:(ld-1)*tres;
ylmin = min(psth.data);
ylmax = max(psth.data);

% Draw the intervals - background rectangles
if isfield(dopts,'intervals')
    intervals = dopts.intervals;
    for i=1:length(intervals)
        hp = patch([intervals(i).ts intervals(i).te intervals(i).te intervals(i).ts],[ylmin ylmin ylmax ylmax],intervals(i).color);
        set(hp,'EdgeColor',intervals(i).color);
    end
end

hold on
if exist('ax','var') & ~isempty(ax)
    hold(ax(2),'on');
end
[oax,hg] = mzplotyy(t,psth.data,psth.contrib,gcolor,contrib_color,ax);
%%scs take out contribution plot

hold off
%hold(oax(2),'off'); %scs take out contrib
xlabel('Time [s]');
if clust >= 0
    ylabel([num2str(elec) '(' num2str(clust) ') ' ylab]);
else
    ylabel([num2str(elec) ' ' ylab]);
end

% Draw any lines specified in dopts and associated text
if isfield(dopts,'markers') & ~isempty(dopts.markers)
    markers = dopts.markers;
    for i=1:length(markers)
        h = line([markers(i).t markers(i).t], [ylmin ylmax]);
        set(h,'Color',markers(i).color);
    end
    text(markers(i).t,ylmax*0.97,markers(i).text);
end

if isfield(dopts,'avbgn')
    avbgn = dopts.avbgn;
    avend = dopts.avend;
    for i=1:length(avbgn)
        nbg = round(avbgn(i) / tres + 1);
        nen = round(avend(i) / tres + 1);
        avp = mean(psth.data(nbg:nen));
        savp = sprintf('%-5.1f',avp);
        text(avbgn(i),avp,savp);
    end
end

axis tight
syl = ylim;
syl(1) = 0;
ylim(syl);
%axes(oax(2)) %scs take out contrib
axis tight
syl = ylim;
syl(1) = 0;
ylim(syl);
%axes(oax(1))%scs take out contrib
