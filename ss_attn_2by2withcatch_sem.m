 clear all
%Attention Task, just a few files --to keep analysis simpler for today

%March 18 2013 - scs
%this is the most up to date analysis-- correct BRT timing
%and can handle different reward values. needs to be updated to accomodate
%catch trials.

%12-4-12 scs see psth_allbrt for all files. 

%PITd
%April 11 2013

load '/Freiwald/ppolosecki/lspace/quincy/130411Quincy/proc/Quincy_130411_0056.mat';
ifname = '/Freiwald/ppolosecki/lspace/quincy/130411Quincy/proc/130411_quincy_011_r.dh5';
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130411Quincy\proc\Quincy_130411_0056.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130411Quincy\proc\130411_quincy_011_r.dh5'
%unit 2-- distractor supression, saccade tuning for 135
%unit 3-- distractor suppression,saccade tuning for 135
% 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130411Quincy\proc\Quincy_130411_0050.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130411Quincy\proc\130411_quincy_006_r.dh5';
%unit 1-- shows suppression of distractor, attention shift (1,2), saccade tuning 135 
%unit 2 -- strong facilitation and saccade tuning (1, 135), 

%April 10 2013
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130410Quincy\proc\Quincy_130410_0055.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130410Quincy\proc\130410_quincy_004_r.dh5';
% 2 units, both good, unit 2 great

%april 9th 2013
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130409Quincy\proc\Quincy_130409_0039.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130409Quincy\proc\130409_quincy_003_r.dh5';
%3 units

%April 5th
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130405Quincy\proc\Quincy_130405_0024.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130405Quincy\proc\130405_quincy_004_r.dh5';

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130405Quincy\proc\Quincy_130405_0027.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130405Quincy\proc\130405_quincy_006_r.dh5';

%April 3, 2013
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130403Quincy\proc\Quincy_130403_0032.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\130403Quincy\proc\130403_Quincy_007_r.dh5'
% unit 1,
%weird, might have attentional inhibtion

% %PITd
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130402Quincy\proc\Quincy_130402_0024.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\130402Quincy\proc\130402_Quincy_006_r.dh5';

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130326Quincy\proc\Quincy_130326_0006.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\130326Quincy\proc\130326_Quincy_001_r.dh5';

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\130326Quincy\proc\Quincy_130326_0017.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\130326Quincy\proc\130326_Quincy_002_r.dh5';
%note can't get surface 3 or 4, just 1 and 2 --FIX scs 3/29/13

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\120508Quincy\proc\120508Quincy_0012.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\120508Quincy\proc\120508_Quincy_001_r.dh5';
%crashes--fix

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\120504Quincy\proc\120504Quincy_0005.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\120504Quincy\proc\120504_Quincy_001_r.dh5';
%very few trials

% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\120504Quincy\proc\120504Quincy_0007.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\120504Quincy\proc\120504_Quincy_002_r.dh5';
%good--used for figures

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\111111Michel\proc\111111Michel_0016.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\111111Michel\proc\111111_michel_001_r.dh5';
% %not so strong attn effect

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\111111Michel\proc\111111Michel_0019.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\111111Michel\proc\111111_michel_002_r.dh5';
%also not so strong


%QUINCY Dec 6 2012 A+3 M+4 (from center, non offset, 0 degree grid, groove
%aligned to anterior screw
% load 'C:\Users\Sara\Desktop\DataAnalysis\Quincy\121206Quincy\proc\Quincy_121206_0035.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Quincy\121206Quincy\proc\121206_quincy_004_r.dh5';
%Strong attention effect
%post saccadic response somewhat consistent, cell shows very little
%saccadic response in MGS, but strong visual response. 


%Nov 27 2012 P + 6 M +2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121127Michel\proc\Michel_121127_0032.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121127Michel\proc\121127_michel_005_r.dh5';
%NOTE, USING 4 SURFACES FOR THIS SESSION
%(Y, unit 1)

%Nov 21 2012 P + 6 M +2

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121121Michel\proc\Michel_121121_0021.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121121Michel\proc\121121_michel_007_r.dh5';
%unit 1, weird attention effect, look more closely-- especially saccade
%effect
%Y unit 2, assuming foveal RF


%Nov 20 2012 P +6 M +2

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\Michel_121120_0023.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\121120_michel_004_r.dh5';
%hico loco, saccade tuning-- analyze further--strong sac tuning for hico,
%stim 1-- more eccentric if attention ref frame. look at remapped response.
%no clear attention effect, might be because surfaces are NOT in RF
%Y!
%Nice example of how postsaccadic (remapping signal?) is also tuned and
%affected, here the surfaces are pretty much outside of the RF
%also very good example
%cue period: opposite attn response if any, converge at point surfaces
%appear

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\Michel_121120_0027.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\121120_michel_006_r.dh5';
%%SPECIAL CASE two target locations, surfaces in RF (eccentric)

%Nov 16 P= 6 M+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\Michel_121116_0028.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\121116_michel_011_r.dh5';
% %Y LIP NEURON-- interesting candidate for attentional coordinate frame
%hypothesis!!!
%Unit 1, example neuron-- attn effect, clear saccade attn effect, use 80%
%coherence
%not motion selective based on tuning analysis
%no attention effect for cue
%hicoloc0 version of task -- interesting, ramping for surf outside rf, attn
%interaction for low coh (attn effect in rf v out is larger than for hico)
%see ramping for 225 targets, some difference in coherence (small)
%******************************************************************
% 
% Feb 15 2013, P + 6 L+ 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130215Michel\proc\Michel_130215_0032.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130215Michel\proc\130215_michel_006_r.dh5';
%note two units, two is better than 1
%unit 1: very weak attn effect if any for surface 1, mild cue effect for
%surface 1
%unit 2 : some cue effect for stim 1, but clear attn effect
%rew: unit 1, not much difference, higher for low rew if anything


% % Feb 14 2013, P + 6 L+ 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130214Michel\proc\Michel_130214_0033.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130214Michel\proc\130214_michel_008_r.dh5';
% %note two units, unit 1 is good
 %unit 1, no cue effect, strange delay (attn) effect, significant early but
 %decays to equal levels late in trial (stim 3 greater than 2)
 %no clear reward modulation, most trials have similar reward though


% % Feb 11 2013, P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130212Michel\proc\Michel_130212_0037.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130212Michel\proc\130212_michel_007_r.dh5';
%unit 1, stim 1 and 3
%marignal effect, too few trials in some conditions to be sure-- not
%reorienting in this case, sac response stronger in possibly attended
%condition
%unit 1-- some attn effect emerges during delay (stronger for 1 than 3), no
%cue effect

% % Feb 11 2013, P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130211Michel\proc\Michel_130211_0030.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130211Michel\proc\130211_michel_010_r.dh5';
%unit 1, cell had no visual response, but showed saccade tuning. No effect.
%(stim 2 and 3)
%unit 1, no cue effect, no attn effect

% Feb 8 2013, P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130208Michel\proc\Michel_130208_0039.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130208Michel\proc\130208_michel_007_r.dh5';
%unit 1, isolation in and out-- but includes neutral trials with surface in
%rf as a gauge of isolation quality
%also note, surfaces are 1 and 3, 2 is the neutral surface
%unit 1, no cue effect, no attn effect

% Feb 6 2013, P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130206Michel\proc\Michel_130206_0020.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130206Michel\proc\130206_michel_005_r.dh5';
%unit 1, isolation in and out.
%unit 1, no cue effect, no attn effect

% % Feb 5 2013, P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130205Michel\proc\Michel_130205_0026.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130205Michel\proc\130205_michel_005_r.dh5';
%unit 1, shows shift in tuning, but also a surface impinges on RF (attn
%effect overall)
%unit 1, shows late emerging attn response in cue period(?) check timing.
%also late attn mod in delay
%USE THIS TO CHECK TIMING OF INTERVALS!! scs mar 3 2013

%Nov 27 2012 P + 6 M +2 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121127Michel\proc\Michel_121127_0032.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121127Michel\proc\121127_michel_005_r.dh5';
%NOTE, USING 4 SURFACES FOR THIS SESSION
%(Y, unit 1)

%Nov 21 2012 P + 6 M +2
% 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121121Michel\proc\Michel_121121_0021.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121121Michel\proc\121121_michel_007_r.dh5';
%unit 1, weird attention effect, look more closely-- especially saccade
%effect
%Y unit 2, assuming foveal RF
%hico loco, not much difference overall (unit 2)hint of seperation for
%180, low firing rate, not many trials so hard to say

%Nov 20 2012 P +6 M +2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\Michel_121120_0023.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\121120_michel_004_r.dh5';
%hico loco, saccade tuning-- analyze further--strong sac tuning for hico,
%stim 1-- more eccentric if attention ref frame. look at remapped response.
%no clear attention effect, might be because surfaces are NOT in RF
%Y!
%Nice example of how postsaccadic (remapping signal?) is also tuned and
%affected, here the surfaces are pretty much outside of the RF
%also very good example
%cue period: opposite attn response if any, converge at point surfaces
%appear

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\Michel_121120_0027.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121120Michel\proc\121120_michel_006_r.dh5';
%%SPECIAL CASE two target locations, surfaces in RF (eccentric)

%Nov 16 P + 6 M + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\Michel_121116_0014.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121116Michel\proc\121116_michel_004_r.dh5';
%unit 1 LIP, strange-- attention effect reversed?? check this--yes, weird
%Y
%reward: actually see higher response for lower reward in presac epoch

% %Nov 15 P+6 M+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121115Michel\proc\Michel_121113_0086.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121115Michel\proc\121115_michel_009_r.dh5';
%unit 1, overall attn effect + ispi, but map is contra and sac effect is
%contra
%Y, but weak, look at map
%hicoloco-- overall, clear ramping for 80 v 30 coherence-- esp aligned to
%brt, difference for 180 target aligned to brt and esp to saccade, does
%show dip at 200ms after stim onset


%Nov 13 2013 P+6 M+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121113Michel\proc\Michel_121113_0046.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121113Michel\proc\121113_michel_007_r.dh5';
%no clear effect cue effect, RF map, faint contra, surface onset contra,
%but attention is ipsi
%hico loco -overall slight difference in coherences (aligned to saccade),
%slight effect for 270,  


%__Nov 8 2012 P+6 L+2

%  load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121108Michel\proc\Michel_121108_0017.mat';
%  ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121108Michel\proc\121108_michel_003_r.dh5';
%two units, unit 1 meager attn, confusing saccade
% (Y)  % unit 2, weak attn effect, saccade tuning

%__ %nov 7 2012 
%  load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121107Michel\proc\Michel_121107_0026.mat';
%  ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121107Michel\proc\121107_michel_005_r.dh5';
 %Y, look at map, unclear)%unit 1, saccade stuff, look closer
% 
%  load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121107Michel\proc\Michel_121107_0041.mat';
%  ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121107Michel\proc\121107_michel_010_r.dh5';
%Y, ipsi effect! %unit 1, meager attn effect overall, some saccade tuning

%__Oct 24 2012
%  load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121024Michel\proc\Michel_121024_0053.mat';
%  ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121024Michel\proc\121024_michel_005_r.dh5';
%Y, but unclear geometry


%__Oct 18 2012 P + 5 L + 3
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121018Michel\proc\Michel_121018_0043.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121018Michel\proc\121018_michel_007_r.dh5';
%no attention effect, saccade weirdness
%Y, but strange geometry--very eccentric lower field?

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121018Michel\proc\Michel_121018_0031.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121018Michel\proc\121018_michel_003_r.dh5';
%meager attentional effect, some saccade tuning
%Y, and good cell isolationwise


%__Oct 17 2012 P + 5 L + 3
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121017Michel\proc\Michel_121017_0035.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121017Michel\proc\121017_michel_008_r.dh5';
%strong attention effect, some saccade tuning, good cell
%Y, good example for attention, RF is upper field
%no cue effect

%__Oct 12 2012 P + 10 L + 1
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121012Michel\proc\Michel_121012_0014.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121012Michel\proc\121012_michel_002_r.dh5';
%saccade tuning, attentionally modulated
%Y consistent with a small, parafoveal RF


%__Oct 4 2012 P + 7 L + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121004Michel\proc\Michel_121004_0037.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121004Michel\proc\121004_michel_001_r.dh5';
%stable, boring cell-- no neutrals condition, high coherence
%Y, maybe if very eccentric upper field

%__Oct 3 2012 P + 7 L + 1
% 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121003Michel\proc\Michel_121003_0050.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121003Michel\proc\121003_michel_006_r.dh5';
%may have some perisaccadic tuning
%Y if lower left field

%__Oct 1 2012 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121001Michel\proc\Michel_121001_0028.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121001Michel\proc\121001_michel_005_r.dh5';
% %no clear effect
%Y-- maybe, look again
%hicoloco overall, some difference in ramping around saccade, less clear
%aligned to brt, rt diffs look modest (check)

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\121001Michel\proc\Michel_121001_0034.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\121001Michel\proc\121001_michel_008_r.dh5';
%LIP neuron, clear tuning
%Y, geometry not clear though


%__Sept 26 2012 P + 7 L + 1
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120926Michel\proc\Michel_120926_0037.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\120926Michel\proc\120926_michel_007_r.dh5';
%LIP neuron, strong attention effect, good saccade effect--nice example
%Y!
%rew- might have slightly higher f.r. for higher rew anticipation, but not
%ramping--constant offset

%Sept 24 2012__P+6 L+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120924Michel\proc\Michel_120924_0025.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\120924Michel\proc\120924_michel_007_r.dh5';
% % two good stable units (1 and 2)
%clear distinction around time of saccade for in/out (trustees poster
%neuron)
%Y unit 1, ipsi upper field, shows attn and saccade effect, but rho was wrong
%for some trials
%Y unit 2, ipsi, clear RF, also good example but rho wrong for some trials,
%no attn effect, but good saccade effect
%rew: unit 1, some effect for higher v lowest, but not clear
%rew: unit 2, no systematic effect

%__Sept 20 2012 P+6 L+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120920Michel\proc\Michel_120920_0006.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\120920Michel\proc\120920_michel_003_r.dh5';
%unit 1, saccade effect, strange-- ipsi?
%M--timing of saccade signal strange, shows effect but not including now.

%__Sept 12 2012 P+6 L+2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120912Michel\proc\Michel_120912_0010.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\120912Michel\proc\120912_michel_001_r.dh5';
%some saccade effect
%Y, good example
%rew--small diff in rew, small diff in neural resp. ns

% 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120906Michel\proc\Michel_120906_0021.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120906Michel\proc\120906_michel_004_r.dh5';
%above file has very good isolation of unit 1, saccadic repsonse and many
%trials, also shows direction tuning
%Y, but can't figure out geometry...
%reward, no clear effect but only two not very different levels
% __

% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120827Michel\proc\Michel_120827_0019.mat'; %UNIT 2 SAC EFFECT
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120827Michel\proc\120827_michel_005_r.dh5';
%Unclear-- RF not capture in MGS, possibly more foveal 135 deg location?


% __
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120823Michel\proc\Michel_120823_0024.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120823Michel\proc\120823_Michel_001_r.dh5';
%cluster 0, 2, attention effect,saccade interaction **,
%(Y?, if RF is foveal)

% 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120725Michel\proc\120725Michel_0061.mat';  %this is the sac/attn interaction
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Michel\120725Michel\proc\120725_Michel_013_r.dh5';
%cluster 0, note to use unsorted spikes (cluster 0), strong attention effect, saccade
%effect **
%Y
%not a clear reward effect, reversed if anything


% TWO BY TWO MICHEL

%April 15 2013 P + 6 L + 2
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130415Michel\proc\Michel_130415_0011.mat';
% ifname = 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130415Michel\proc\130415_Michel_005_r.dh5';
%unit 1, very few trials, attention effect
%
%April 16 2013 P + 6 L + 2  
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130416Michel\proc\Michel_130416_0030.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130416Michel\proc\130416_michel_007_r.dh5';
  %unit 1, small attn effect--note, not many trials, no clear saccade
  %tuning (but not in MGS either)
  %* unit 2, --strong attention effect, not many trials
  %2 UNITS

%April 17 2013 p + 6 L + 2
 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130417Michel\proc\Michel_130417_0021.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130417Michel\proc\130417_michel_006_r.dh5';
%unit 1, strong attention effect--covert throughout (compare 135 brt for surf 1 v 3 v 2&4), tuned for all saccades (same as mgs, no clear ramping) -- decent number of trials
%unit 2, ipsi? attn effect, tuned for all saccades, a little weird.   

%April 18 2013 p + 6 L + 2
 
% load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130418Michel\proc\Michel_130418_0017.mat';
% ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130418Michel\proc\130418_michel_004_r.dh5';
%unit 1, some attention effect, upper field RF
%unit 2, some attetnion effect, weird. Both cells had poor isolation.   
 
%April 19 2013 p + 6 L + 2
 
%load 'C:\Users\Sara\Desktop\DataAnalysis\Michel\130419Michel\proc\Michel_130419_0015.mat';
%ifname ='C:\Users\Sara\Desktop\DataAnalysis\Michel\130419Michel\proc\130419_michel_007_r.dh5';

%%

addpath(fullfile('..','akinori_analysis'))
addpath(fullfile('..','psth_indiv_trials'))
addpath(fullfile('..','detlef_psth'))
addpath('/rdata/space/opt/local_matlab_tools/dhfun3/')
c_list= [[1 1 1]; [0 1 1]; [0 0 1];[0 0 0]; [0 1 0]; [1 1 0]; [1 0 0]; [1 0 1]];
bid = 1; % Channel number (1 for Blackrock recordings in Freiwald lab)
tres = 1e9/30000; % Specified in nanoseconds, reciprocal of sample rate

seloutcomes = [7]; % Trial outcomes to analyze (Visiko outcome 7 = success)

clussel =2; % Cluster number(s)

%Attention task choices
selstimnos = [1];%1; % One stimulus number per graph (attended 2/nonattended condition 1)
selbrt = []; %0,45 180, 225 %select the BRT for heiko task --which directions do you want to see? scs
selcoh = [];% any value-- precise or not-- prompt will appear if not empty
%80; %25 30 and 80 are the choices for now 
rewval = [];

allbrt =0; % this is to automate the process for all brt's %fix to collapse

align_brt =0; % to do this, you must also select to align to the fixation out.  

align = [2]%1 is align to stimulus onset, 2 is align to the fixation out, 3 is cue only


overlay =0;%for the first plot, be sure to set to 0
overlaycolor = 'r'; 
type = 1; % 1 is normal, 2 is neutral

SEM = 1; %adds SEM to psth's, also adds a lot of time. 

% to look at neutral trials, all selection criteria must be empty (except for stimnos), allbrt
% 0, align = 1, type = 2
%NOTE that surf_str doesn't grab info for neutral trials--fix SCS feb 12
%2013

%[task_fig] = task_config(ifname, tds) 

%___________________________________________________________________________________________
if align == 1
start_trig.type ='Marker'%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig.name ='Stimulus_onset'%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig.offset = .8%-.05%1.0%-.25%-.05; %usually a negative value 
end_trig.type ='Marker'%'TrialEnd'; %'Marker'%%
end_trig.name ='Stimulus_onset'%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig.offset =2.5%1.8;% 4.0%.05;%3.0
end

%just cue period

if align == 3
start_trig.type ='Marker'%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig.name ='Stimulus_onset'%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig.offset = -.05%-.05%1.0%-.25%-.05; %usually a negative value 
end_trig.type ='Marker'%'TrialEnd'; %'Marker'%%
end_trig.name ='Stimulus_onset'%'Fixation_out'%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig.offset =1;%1.8;% 4.0%.05;%3.0
end

%full trial
% if align == 1
% start_trig.type ='TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
% start_trig.name ='TrialStart';%;%'Stimulus_onset';
% start_trig.offset = -.05%-.05; %usually a negative value 
% end_trig.type ='TrialEnd'; %'Marker'%%
% end_trig.name ='TrialEnd';
% end_trig.offset = .05%.05;
% end

if align ==2
start_trig.type ='Marker'%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig.name ='Fixation_out'%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig.offset =-1.0%-1.0%-2.0; %usually a negative value (offset doesn't work when you align brt, fix)
end_trig.type ='Marker'%'TrialEnd'%'TrialEnd'; %'Marker'%%
end_trig.name ='Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig.offset = 0.3%.15;
end

psthfilter.type = 'Gauss';
psthfilter.fs        = round(1e9/tres); %Sampling frequency (Hz)
psthfilter.sigma_ms  = 25; %Specified in milliseconds
psthfilter.length_ms = 100;   %Specified in milliseconds, truncation length, usually 4x more than sigma
% start_trig.type ='Marker'%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
% start_trig.name ='Stimulus_onset'%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
% start_trig.offset = -1.0; %usuallty a negative value (offset doesn't work when you align brt, fix)
% end_trig.type ='Marker'%'TrialEnd'; %'Marker'%%
% end_trig.name ='Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
% end_trig.offset = .5;
psthdec = 100;

% ------------ END OF USER-EDITABLE BLOCK -----------------
clear dhfun
% Design the filter
psthfilter = design_psth_filter(psthfilter);


fid = dhfun('open',ifname);
tmap = dh_get_trialmap_struct(fid);
trls = [tds{1}.trials]; %holds all the basic trial parameters

%this gets you the name, number, rho and phi for the target surface of
%every trial
%[surf_str] = surface_info_withcatch(tmap,tds);
[surf_str] = surface_info_withcatch(tds);
surface_numbers = [surf_str.numb];


%get the info for the heiko task
    td = tds{1};
    nocue = false;
    p = heiko_log_surf_params(td,nocue);
    p.flickerdur;
    p.pausedur;


%this is for the special case of having more than two surfaces, otherwise
%tmap does not get changed. 
%if length(unique(surface_numbers))==2 %this is stupid, fix it SCS feb 6 2013
    tmap.sn = [surf_str.numb]';
%end

% Apply trial selection
for j= 1:length(tds{1}.trials)
    if ~strcmp(tds{1}.trials(j).type,'Normal') brt(j) = nan; 
    else
    brt(j)= tds{1}.trials(j).targsurf.brtdir;
    end 
end
brtz = unique(brt); 
brts = (brtz(~isnan(brtz)));

% for j= 1:length(tds{1}.trials)
%     if type ==1; 
%         if strcmp(tds{1}.trials(j).type,'Neutral')
%         end
%         if type == 2 ; 
%         end
        


if allbrt ==0;
    %this runs the code just once for all the brt's collapsed, but it is
    %the same code that gets used below for the loop. find a way to make
    %more effecient. scs
 
    %tmap = dh_get_trialmap_struct(fid); % restore tmap for each iteration
if ~isempty(selbrt)
    tmap= trialselect_brt(tmap, selbrt, tds); % this must go first, it can't handle less than the full trial map--fix this! see coh scs
end
if ~isempty(selcoh)


all_cohs = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.brtCoher];
coh_vals = str2double(unique(all_cohs));
if isempty(find(selcoh == coh_vals)) && ~isempty(selcoh)
     coh_vals
    selcoh = input('Choose a valid coherence value from above  ') 
end
tmap = trialselect_coh(tmap, selcoh, tds); %select trials based on coherence level of BRT
end
% if ~isempty(seltargets) %select the targets of choice% added this to grab specific targets scs
% tmap = trialselect_targets(tmap, seltargets,t_idx); 
% end

if nansum(tmap.sn)>0 && ~isempty(selstimnos);%scs %select the stimulus (surface)
tmap = trialselect_stimnos(tmap,selstimnos);
end
if nansum(tmap.oc==seloutcomes) >0; %scs %select outcomes of trial (success)
tmap = trialselect_outcomes(tmap,seloutcomes);
end

%select here for reward value, make this work and give values to choose
%from like selecting the coherence. scs march 12 2013
if ~isempty(rewval)
    rvs = [];
    rvs = [surf_str.rew];
    allrews = unique(rvs);
if ~isempty(rewval) && isempty(find(rewval == allrews))
    allrews
    rewval = input('Choose a valid reward value from above   ')
end
    tmap = trialselect_rewvals(surf_str,rewval,tmap);
end

    %new bit to make sure catch trials do not sneak in
    if type ==1
        for j = 1:length(trls);
            if strcmp(tds{1}.trials(j).type, 'Normal') ==1
                normal(j) = 1;
            else
                normal(j) = 0;
            end
        end
    end
    
        tmap = trialselect_normal(normal, tmap); % include in function
    

%this gives you coherences of BRT
%so if you want to align to the brt, add the appropriate amount of offset
%time to the interval of interest

% if align_brt ==1;
%     %note, here is code to do what you already did below...but standardized
% 
%     for n = 1:length(tmap.idx)
%        Brt_time(n) = tmap.ts(n) + (p.pausedur(tmap.idx(n))*1000000)+(p.bit_total_time(tmap.idx(n))*1000000);
% %         brt_time(n) = tmap.ts(n) + (p.pausedur(tmap.idx(n))*tres)+(p.bit_total_time(tmap.idx(n))*tres);
%     end
% end
%surf_angle = [p.cuedsurfangles(1,:)]; %this is the phi for the surface that was cued for each trial


% Apply the reference points %scs note that this redefines ts and te to be
% intervals start and end, not trial start and end-- verify this. 
[ts,te,inclidx] = find_ref_points(fid,tmap,start_trig,end_trig);
%new bits scs 12-10-15
% inclidx = inclidx(ts>0);
% ts = ts(ts>0);
% te = te(ts>0); 
%if your ts is negative, everything crashes. eliminate negative ts's here!
if find(ts<0)
    disp('Hey, your ts is less than zero in this category.') 
end




%%
% % try to get BRT time independent of triggers selected

start_trig2.type ='Marker';%;%'TrialStart' %'Marker';%%'Marker';%'TrialStart'; %'Marker';
start_trig2.name ='Stimulus_onset';%'Fixation_out';%'TrialStart';%;%'Stimulus_onset';
start_trig2.offset = 0;%-.05%1.0%-.25%-.05; %usually a negative value 
end_trig2.type ='TrialEnd'; %'Marker'%%
end_trig2.name ='TrialEnd';%'Stimulus_onset'%'Fixation_out'%'TrialEnd'; %'Fixation_out'%'TrialEnd';
end_trig2.offset =0;%1.8;% 4.0%.05;%3.0


 [Stimon_time,Trlend_time,inclidx2] = find_ref_points(fid,tmap,start_trig2,end_trig2);
    pause = p.pausedur(tmap.idx(1))*1e9;
   stimon_time = Stimon_time(inclidx);
   trlend_time = Trlend_time(inclidx); 
    for n = 1:length(tmap.idx(inclidx))
        brt_time(n) = stimon_time(n) + pause+((p.bit_total_time(tmap.idx(inclidx(n))))*1e6);
%         brt_time(n) = tmap.ts(n) + (p.pausedur(tmap.idx(n))*tres)+(p.bit_total_time(tmap.idx(n))*tres);
    end

    
    rts = [];
    Rts = [];
    Rts = [te-brt_time'];
    rts = Rts/1e6;
    mean_rt = nanmean(rts);
    std_rt = nanstd(rts);
    num_rt = length(rts);
    sem_rt = std_rt/sqrt(num_rt); 
%%

    
    
    
%quick substitution here scs
if align_brt ==1
    ts = brt_time';
end


% Load the spikes
nspk = double(dhfun('GETNUMBERSPIKES',fid,bid)); %the numbers of spikes (total)
tstamps = dhfun('READSPIKEINDEX',fid,bid,1,nspk); %times that every spike occured

[spike_clus] = dhfun('READSPIKECLUSTER',fid,bid,1,nspk); %tells you cluster information
spks = tstamps(spike_clus == clussel); 
 



% Calc the PSTH
clear psth
if length(tmap.idx)>0
psth = calc_psth(spks,ts,te,ts,tres,psthfilter.Num,1,psthfilter.delay,psthdec);
if SEM  ==1; 
[psth_se] = dh_calc_psth_stderr(fid,bid,ts,te,clussel,tres,psthfilter.Num,psthfilter.Den,psthfilter.delay,psthdec);
sem = [];
sem = psth_se.stderr;
mean = [];
mean = psth.data;
mean_idx = 1:length(mean); 
upperlim = [];
lowerlim = [];
upperlim = mean + sem;
lowerlim = mean - sem;

Tres = tres / 1e9; % Make tres be in seconds
ld = length(psth.data);
T = 0:Tres:(ld-1)*Tres;
t = T *psthdec;
end
    % And plot it
    if overlay ==1
        hold on
        if SEM ==1
            legendlist = ['test']
            color2 = overlaycolor; % [1 0 0];
            [maxvalue minvalue]= shadowcaster_ver3(t', mean', sem', legendlist,color2) %psth_plot(psth,tres*psthdec,[],[1 0 0],[1 0 1],[],bid,clussel,'Spikes per second'); %can add intervals to empty field
    hold on
        end
        psth_plot(psth,tres *psthdec,[],overlaycolor,[1 0 1],[],bid,clussel,'Spikes per second'); %can add intervals to empty field
    else
        figure
        %whitebg([0.9 0.9 0.9]);
        if SEM ==1
             legendlist = ['test']
             color1 = [0 0 1]
            [maxvalue minvalue]= shadowcaster_ver3(t', mean', sem', legendlist,color1) %psth_plot(psth,tres*psthdec,[],[1 0 0],[1 0 1],[],bid,clussel,'Spikes per second'); %can add intervals to empty field
        hold on
        end
        psth_plot(psth,tres *psthdec,[],[0 0 1],[1 0 1],[],bid,clussel,'Spikes per second'); %note that you can get a figure for each brt
    end
 
    
    % scs Make the time axis relative to the onset of the stimulus
    % axis[(start_trig.offset )]
    hold on
  line([0-start_trig.offset 0-start_trig.offset], ylim,'Color', 'k', 'LineStyle', '-.') 
  
  if align ==1
  line([(0-(start_trig.offset) +p.pausedur(1)) (0-start_trig.offset +p.pausedur(1))], ylim, 'Color', 'r', 'LineStyle', '--')
 hold on
 %xlim([0.05 2])
  end
  
    if ~isempty(selbrt)
    Brts = int2str(selbrt);
    legend(Brts,'Location', 'Northwest')
    end
end
end


%this iterates the code for each brt
if allbrt ==1;
    figure
   %whitebg([0.9 0.9 0.9]);
    overlay =1; % comment this is if you want a figure for each BRT
%select the brt of interest for Heiko task
for b = 1:length(brts)
    selbrt = brts(b);
    tmap = dh_get_trialmap_struct(fid); % restore tmap for each iteration

   
      
if ~isempty(selbrt)
    tmap= trialselect_brt(tmap, selbrt, tds); % this must go first, it can't handle less than the full trial map--fix this! see coh scs
end   
    %this is for the special case of having more than two surfaces, otherwise
%tmap does not get changed. 
if length(unique(surface_numbers(find(~isnan(surface_numbers)))))==2 %double check this SCS Feb 6 2013
    tmap.sn = [surf_str.numb]';
end
    



if ~isempty(selcoh)
    
all_cohs = [tds{1}.doc_data(1).FIELDPARAMS_HEIKO.brtCoher];
coh_vals = str2double(unique(all_cohs));
if isempty(find(selcoh == coh_vals)) && ~isempty(selcoh)
     coh_vals
    selcoh = input('Choose a valid coherence value from above  ') 
end
tmap = trialselect_coh(tmap, selcoh, tds); %select trials based on coherence level of BRT
end

%grab trials with a specific reward value
if ~isempty(rewval)
    rvs = [];
    rvs = [surf_str.rew];
    allrews = unique(rvs);
if ~isempty(rewval) && isempty(find(rewval == allrews))
    allrews
    rewval = input('Choose a valid reward value from above   ')
end
    tmap = trialselect_rewvals(surf_str,rewval,tmap);
end

if nansum(tmap.sn)>0 && ~isempty(selstimnos);%scs %select the stimulus (surface)
tmap = trialselect_stimnos(tmap,selstimnos);
end
if nansum(tmap.oc==seloutcomes) >0; %scs %select outcomes of trial (success)
tmap = trialselect_outcomes(tmap,seloutcomes);
end


%this gives you coherences of BRT
%so if you want to align to the brt, add the appropriate amount of offset
%time to the interval of interest

if align_brt ==1;
    %note, here is code to do what you already did below...but standardized
    td = tds{1};
    nocue = false;
    p = heiko_log_surf_params(td,nocue);
    p.flickerdur;
    p.pausedur;
    %this might need to be fixed
    % pause = p.pausedur(tmap.idx(1))*1e9
    
    for n = 1:length(tmap.idx)
        brt_time(n) = tmap.ts(n) + (p.pausedur(tmap.idx(n))*1000000)+(p.bit_total_time(tmap.idx(n))*1000000);
    end
end
% Apply the reference points %scs note that this redefines ts and te to be
% intervals start and end, not trial start and end-- verify this. 
[ts,te,inclidx] = find_ref_points(fid,tmap,start_trig,end_trig);

%quick substitution here scs
if align_brt ==1
    ts = brt_time(inclidx)';
end


    

% Load the spikes
nspk = double(dhfun('GETNUMBERSPIKES',fid,bid)); %the numbers of spikes (total)
tstamps = dhfun('READSPIKEINDEX',fid,bid,1,nspk); %times that every spike occured

[spike_clus] = dhfun('READSPIKECLUSTER',fid,bid,1,nspk); %tells you cluster information
spks = tstamps(spike_clus == clussel); 


% Calc the PSTH
clear psth
if tmap.idx>0
    psth = calc_psth(spks,ts,te,ts,tres,psthfilter.Num,1,psthfilter.delay,psthdec);
    if SEM ==1
    [psth_se] = dh_calc_psth_stderr(fid,bid,ts,te,clussel,tres,psthfilter.Num,psthfilter.Den,psthfilter.delay,psthdec);
sem = [];
sem = psth_se.stderr;
mean = [];
mean = psth.data;
mean_idx = 1:length(mean); 
upperlim = [];
lowerlim = [];
upperlim = mean + sem;
lowerlim = mean - sem;
Tres = tres / 1e9; % Make tres be in seconds
ld = length(psth.data);
T = 0:Tres:(ld-1)*Tres;
t = T *psthdec;
    end
    % And plot it
    if overlay ==1
        if SEM ==1
        colorlist = [c_list(b,:)];
        lengendlist = ['test'];
        [maxvalue minvalue]= shadowcaster_ver3(t', mean', sem', lengendlist,colorlist);
        hold on
        end
        hold on
        psth_plot(psth,tres *psthdec,[],c_list(b,:),[1 0 0],[],bid,clussel,'Spikes per second'); %can add intervals to empty field
    else
        figure
       %whitebg([0.9 0.9 0.9]);
        if SEM ==1
          colorlist = [c_list(b,:)];
        lengendlist = ['test'];
        [maxvalue minvalue]= shadowcaster_ver3(t', mean', sem', lengendlist,colorlist);
        hold on
        end
        psth_plot(psth,tres *psthdec,[],c_list(b,:),[0 0 1],[],bid,clussel,'Spikes per second'); %note that you can get a figure for each brt 
    end
  %(psth,tres,dopts,gcolor,contrib_color,ax,elec,clust,ylab)

% scs Make the time axis relative to the onset of the stimulus
% axis[(start_trig.offset )]
hold on
line([0-start_trig.offset 0-start_trig.offset], ylim,'Color', 'k', 'LineStyle', '-.') 
if align ==1
  line([(0-(start_trig.offset) +p.pausedur(1)) (0-start_trig.offset +p.pausedur(1))], ylim, 'Color', 'r', 'LineStyle', '--')
end

else
    disp('Boink! No trials fit the criteria you selected')
  end
% Brts = int2str(selbrt);
% legend(Brts,'Location', 'Northwest')
end

end

% %figure
% p_rad = deg2rad(brts);
% rhoz(1:length(p_rad)) = deal((tds{1}.doc_data(1).TARGETSPARAMS.posRho));
% rho = str2double(rhoz);

% figure
% for p = 1:length(p_rad)
%    plot(brts(p), rho(p),'o', 'Color', c_list(p,:))
%     hold on
% end
% figure
% for p = 1:length(p_rad)
% polar(p_rad(p),rho(p),'o'),% 'Color',c_list(p,:))
% hold on
% end

%Using the log file to find the BRT time, verify using above method. scs
% for j= 1:length(tds{1}.trials)
%     if strcmp(tds{1}.trials(j).type,'Neutral') brt_dir(j) = nan; 
%     else
%     brt_dir(j)= tds{1}.trials(j).targsurf.brtdir;
%     bit_num(j)=  nansum(tds{1}.trials(j).targsurf.bitdurs);
%     fr = tds{1}.framerate;
%     bit_totaltime_ms(j) = (bit_num(j)*1000)/fr; %check this make the units correct
%     end 
% end
% 
% for n = 1:length(tmap.idx)
%     brt_time(n) = tmap.ts(n) + (bit_totaltime_ms(tmap.idx(n))*1000000);
% end






dhfun('close',fid);




%sanity check-- plots the trials in time, brt should occur between trials
%start and trial end.
% figure
% for k= 1:length(tmap.ts)
% hold on
% plot(k,tmap.ts(k),'ko')
% hold on
% plot(k,tmap.te(k),'bo')
% hold on
% plot(k,brt_time(k),'rx')
% hold on
% plot(k,brt_time_alt(k),'go')
% end