load ET_data 
% The ocular data is assuemed to be converted into .mat file with the name: ET_data 
% It is also assumed to containt the following data arrays:
% EH_gaze_hcoord: Horizontal Gaze Coordinate over time
% EH_gaze_vcoord: Vertical Gaze Coordinate over time
% EH_gaze_length: Perpendicular distance between the calibration surface
% and the participant/user over time
% PplDiam: Pupil diameter over time
% zero and nan mean no detection of the pupil or closed eyes
% cr_diam/cr_det: Corneal diameter over time
% zero and nan mean no detection of the cornea or closed eyes
% scen_num: scene number over time (optional: to detect whether the
% gaze data is refering to the calibrated surface or not
% It is zero for the default scene defined in ASL Eye-trac7

% % optional interpolation for missing values
% EH_gaze_hcoord = fillmissing(EH_gaze_hcoord,'linear');
% EH_gaze_vcoord = fillmissing(EH_gaze_vcoord,'linear');
% EH_gaze_length = fillmissing(EH_gaze_length,'linear');

% % Optional one-sample spike removal
% % Larsson, Linnéa, Marcus Nyström, and Martin Stridh. "Detection of 
% % saccades and postsaccadic oscillations in the presence of smooth 
% % pursuit." IEEE Transactions on biomedical engineering 60, no. 9 (2013): 2484-2493.
% EH_gaze_hcoord = medfilt1(EH_gaze_hcoord,3);
% EH_gaze_vcoord = medfilt1(EH_gaze_vcoord,3);
% EH_gaze_length = medfilt1(EH_gaze_length,3);

Fs = 360; % Sampling frequency (e.g. 360 Hz)
model_pd = 77; % scaling factor for pupil diameter measurement which 
% should be measured seperately
PplDiam = PplDiam.*3.96/model_pd; % Scaling of pupil size from px to mm in ASL Eye-trac7
u_idx = find(PplDiam~=0 & isnan(PplDiam)==0);
b_samples = find(PplDiam==0 | isnan(PplDiam));
closd_smpl = b_samples; % closed eyes samples
D = diff([0,diff(closd_smpl')==1,0]);
ce_on = closd_smpl(D>0);
ce_off = closd_smpl(D<0);
ce_d = (1 + find(D<0) - find(D>0))/Fs;
ce_on(ce_d<0.1 | ce_d>1)=[];ce_off(ce_d<0.1 | ce_d>1)=[]; 
% Acceptable closed eyes duration >100 ms and <1s
closd_smpl = arrayfun(@colon, ce_on, ce_off, 'Uniform', false);
closeys_idx = cell2mat(closd_smpl');

% ---------- Pupil size preprocessing--------------------------------------
win=10;mis_idx=[];
i=win+1;
while i<length(PplDiam)
    k=1;
    while PplDiam(i+k)==0 && (i+k)<length(PplDiam)
        k=k+1;
    end
    if k>1 && i+k+win<length(PplDiam)
        mis_idx = [mis_idx i-win:i+k-1 i+k:i+k+win];
        i=i+win;
    end
    i=i+k;
end
mis_idx=unique([mis_idx find(PplDiam<(median(PplDiam(u_idx))-3*std(PplDiam(u_idx))))']);
PplDiam(mis_idx)=nan;
PplDiam=fillmissing(PplDiam,'linear');
P_fltr = designfilt('lowpassiir', 'FilterOrder', 3, 'HalfPowerFrequency', ...
    2, 'SampleRate', Fs, 'DesignMethod', 'butter');
% Privitera, Claudio M., et al. "Pupil dilation during visual target detection."
% Journal of Vision 10.10 (2010): 3-3.
% cut-off frequency 4 Hz
PplDiam = filtfilt(P_fltr,PplDiam);
% PplDiam((PplDiam>median(PplDiam)+1) |
% (PplDiam<median(PplDiam)-1))=median(PplDiam);
% --------- End of pupil size preprocessing--------------------------------

% --------- Blink detection------------------------------------------------
L = length(EH_gaze_hcoord);
D = diff([0,diff(b_samples')==1,0]);
b_on = b_samples(D>0);
b_off = b_samples(D<0);
b_d = (1 + find(D<0) - find(D>0))/Fs;
b_outrng = b_d<0.1 | b_d>.4; % Blink duration criteria between 100 and 400 ms
% Eye blink can be even shorter in duration (e.g. 50 ms)
% Stern, John A., Larry C. Walrath, and Robert Goldstein. "The endogenous 
% eyeblink." Psychophysiology 21, no. 1 (1984): 22-33.
undef_idx = arrayfun(@colon, b_on, b_off, 'Uniform', false);
undef_idx = cell2mat(undef_idx');
b_on(b_outrng)=[];b_off(b_outrng)=[];b_d(b_outrng)=[];
undef_idx = unique([undef_idx find(scen_num~=0)']);
b_samples = arrayfun(@colon, b_on, b_off, 'Uniform', false);
b_samples = cell2mat(b_samples'); % blinking samples
% -------- End of blink detection------------------------------------------

% -------- saccde detection------------------------------------------------
[~,g] = sgolay(2,19);
dx = zeros(L,2);
dy = zeros(L,2);
dt = 1/Fs;
for p = 1:2
    dx(:,p) = conv(EH_gaze_hcoord, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    dy(:,p) = conv(EH_gaze_vcoord, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end
for i=1:L-1
    Phi(i) = 2*atand(1/(2*(median(EH_gaze_length(i:i+1)))));
end
for i=1:L-1
    dtheta(i) = Phi(i)*sqrt(dx(i,1).^2+dy(i,1).^2); % Gaze Velocity
    ddtheta(i) = Phi(i)*sqrt(dx(i,2).^2+dy(i,2).^2); % Gaze Acceleration
end

dtheta(L)=dtheta(L-1);
ddtheta(L)=ddtheta(L-1);
% Unacceptable/erroneous samples based on velocity and acceleration limits
err_smpl = find(dtheta>500 | ddtheta>50000); 
int_err_smpl=[];
for j=1:length(err_smpl)-1
    if err_smpl(j+1)-err_smpl(j)<20
        int_err_smpl = [int_err_smpl err_smpl(j):err_smpl(j+1)];
    end
end
err_smpl = sort([err_smpl int_err_smpl]);
err_smpl(diff(err_smpl)==0)=[];
D = diff([0,diff(err_smpl)==1,0]);
err_smpl_on = err_smpl(D>0);
err_smpl_off = err_smpl(D<0);
err_smpl_d = (1 + find(D<0) - find(D>0))/Fs;
long_err_smpl=[];
err_smpl_on(err_smpl_d<0.020)=[];err_smpl_off(err_smpl_d<0.020)=[];
long_err_smpl = arrayfun(@colon, err_smpl_on, err_smpl_off, 'Uniform', false);
long_err_smpl = cell2mat(long_err_smpl);

% Initialization of saccade detection algorithm based on estimated baseline
% noise in saccade velocity
% Based on
% Nyström, M. and Holmqvist, K., 2010. An adaptive algorithm for fixation, 
% saccade, and glissade detection in eyetracking data. Behavior research methods, 42(1), pp.188-204.

PT(1) = 100; % initial value for saccade onset detection
PT(2) = mean(dtheta(dtheta<PT(1)))+6*std(dtheta(dtheta<PT(1)));
i=3;
while abs(PT(i-1)-PT(i-2))>1
    PT(i)= mean(dtheta(dtheta<PT(i-1)))+6*std(dtheta(dtheta<PT(i-1)));
    i=i+1;
end
s_onset_th = min(mean(dtheta(dtheta<PT(end))) + 3*std(dtheta(dtheta<PT(end))),45);
%                 s_smpl = find(dtheta>s_onset_th);
s_smpl = find(dtheta>max(s_onset_th,25));
if isempty(s_smpl)==0
    s_smpl(ismembc(s_smpl,b_samples))=[];
end
if isempty(err_smpl)==0
    s_smpl(ismembc(s_smpl,err_smpl))=[];
end
s_smpl(ismembc(s_smpl,find(isnan(cr_diam) | cr_diam==0)))=[];
%     D = [];
D = diff([0,diff(s_smpl)==1,0]);
s_on = s_smpl(D>0); % saccade onset
s_off = s_smpl(D<0); % saccade offset
s_d = (1 + find(D<0) - find(D>0))/Fs; % saccade duration
s_rmv = find(s_d<0.020 | s_d>0.200); % saccade duration criteria
s_on(s_rmv)=[];s_off(s_rmv)=[];s_d(s_rmv)=[];
s_on=s_on';s_off=s_off';
s_pv = arrayfun(@(s, e) max(dtheta(s:e)), s_on, s_off); % saccade peak velocity
s_a = arrayfun(@(s,e) (mean(dtheta(s:e))).*length(s:e)/Fs, s_on, s_off); % saccade amplitude
rmv_idx = find(sqrt(s_pv-(24+26*s_a))>5 | s_a>20); 
% (Optional) saccade removal criteria based on peak velocity-amplitude relationship and saccade amplitude
% Behrens, F., MacKeben, M. and Schröder-Preikschat, W., 2010. An improved
% algorithm for automatic detection of saccades in eye movement data and 
% for calculating saccade parameters. Behavior research methods, 42(3), pp.701-708.
s_pv(rmv_idx)=[];s_a(rmv_idx)=[];s_on(rmv_idx)=[];s_off(rmv_idx)=[];s_d(rmv_idx)=[];
% -------------- End of saccde detection 
undef_idx = unique([undef_idx rmv_idx' s_rmv]);  % index of undefined samples


% --------- Fixation detection---------------------------------------------
f_samples = 1:L;
cr_rmv = find(isnan(cr_diam) | cr_diam==0);
f_rmv = unique([s_smpl,b_samples,long_err_smpl,undef_idx,cr_rmv']);
f_samples(f_rmv)=[];

D = diff([0,diff(f_samples)==1,0]);
f_on = f_samples(D>0);
f_off = f_samples(D<0);
if isempty(f_on)==0
    f_rmv = find((f_on(2:end)-f_off(1:end-1))<5); % Minimum time between adjacent fixations to combine (11 ms in 360 Hz) as the same as minimum saccade duration
    % But this says (75 ms) which is not accurate enough Ref: Olsen, A. (2012). The Tobii I-VT fixation filter. Tobii Technology.
    f1= f_on(1);f2=f_off(end);
    f_on(1)=[];f_off(end)=[];
    f_on(f_rmv)=[];f_off(f_rmv)=[];
    f_on = [f1 f_on];f_off = [f_off f2];
    
    f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
    f_samples = cell2mat(f_samples);
    
    D = diff([0,diff(f_samples)==1,0]);
    f_d = (1 + find(D<0) - find(D>0))/Fs;
    f_rmv2 = find(f_d<0.04 | f_d>2.5); % fixation duration criteria
    f_on(f_rmv2)=[];f_off(f_rmv2)=[];f_d(f_rmv2)=[];
    f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
    f_samples = cell2mat(f_samples);
    
    
    rmv = [];
    for j=1:numel(f_on)
        C_GP_dist = sqrt((EH_gaze_hcoord(f_on(j):f_off(j))-median(EH_gaze_hcoord(f_on(j):f_off(j)))).^2+(EH_gaze_vcoord(f_on(j):f_off(j))-median(EH_gaze_vcoord(f_on(j):f_off(j)))).^2);
        if sum(C_GP_dist>1)>11 % if the distance to the center of the fixation cluster is higher than 1 degree for more than 11 samples
            rmv = [rmv j];
        end
    end
    f_on(rmv)=[];f_off(rmv)=[];f_d(rmv)=[];
    f_samples=[];
    f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
    f_samples = cell2mat(f_samples);
    undef_idx = unique([undef_idx cr_rmv' f_rmv f_rmv2]);
else
    f_d=nan;
end
% -------------- End of fixation detection---------------------------------

% Further considerations:

% Note that some of the criteria and thresholds may be modified
% depending on the task and settings of the recordings.

% The interpolations can be used for visualization but not for the
% computation of oculometrics

% This algorithm parse the gaze data into saccades, fixations, blinks, and
% undefined events. 

% It also provides pre-processings for pupil size

% Post-saccade oscillations are merged into fixations if they satisfy
% fixation criteria

% The detection of the ocular events is a sequential meaning that the
% fixation detection is dependent ot saccade detection