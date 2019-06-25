% ICA: Index of Cognitive Load calculation
% Based on the patent titled:
% Method and apparatus for eye tracking and monitoring pupil dilation to evaluate cognitive activity
% ------------------------------
Fs = 360; % sampling frequency
EM_SF = 77; % Measured scaling factor from Eye Model in Eyetrac7
PplDiam = PplDiam.*3.96/EM_SF; % Pupil size in mm in Eyetrac7
b_idx = [];
b_idx = find(ppl_rec==0 | isnan(ppl_rec)); % the samples where pupil was closed
D = diff([0,diff(b_idx')==1,0]);
b_on = b_idx(D>0);
b_off = b_idx(D<0);
b_d = (1 + find(D<0) - find(D>0))/Fs;
rmvidx=[];outidx=[];rmv=[];
% preprocessing
% Remove pre-blink and post-blink artifacts (70-100 ms on either side of a blink)
% removing any observation that differs by more than 0.1 mm from the
% baseline
for i=1:length(b_on)
    if b_on(i)-38>0
        preb=find(abs(PplDiam(b_on(i)-38)-PplDiam(b_on(i)-37:b_on(i)-1))>.1); % pre-blink samples
    end
    if b_off(i)+38<L
        prob=find(abs(PplDiam(b_off(i)+38)-PplDiam(b_off(i)+1:b_off(i)+37))>.1); % post-blink samples
    end
    rmvidx = [b_on(i)-(38-preb)' b_off(i)+prob' rmvidx];
end
outidx = find(abs(diff(PplDiam))>0.1); % diff between two consecuitive observations
difP = PplDiam(3:end) - PplDiam(1:end-2); % diff between two observations in a three consecutive observations
outidx1 = find(difP>0.1)+2;
% concatenation of invalid samples
rmv = unique([outidx' outidx1']);
rmv = unique([rmv rmvidx]);
rmv = unique([b_idx' rmv]);
rmv(rmv<1)=[];
PplDiam(rmv)=nan;
PplDiam = fillmissing(PplDiam,'linear'); % linear interpolation of invalid samples

% wavelet analysis
% wavelet decomposition
[c,l] = wavedec(PplDiam,3,'db22'); % db22 is chosen to cover approx. 22*2*1000/360=122 ms
XD = wden(c,l,'minimaxi','h','mln',3,'db22'); % denoising using minimax estimation
[c,l] = wavedec(XD,3,'db22'); % decomposing the denoised signal
X1 = wrcoef('d',c,l,'db22',1);
X2 = wrcoef('d',c,l,'db22',2);
X3 = wrcoef('d',c,l,'db22',3);
X=X1+X2+X3; % signal reconstruction


% The following is my suggestion to read a collection of fluctuations as
% one, otherwise this step is not mentioned in the patent
[yupper,~] = envelope(X,100,'analytic');
yy = smooth(yupper);
[pks,locs] = findpeaks(yy,Fs,'MinPeakDistance',0.2,'MinPeakHeight',0.01);
impuls_numbers = numel(locs); % number of pupillary peaks (phasic response)

% One alternative to ICA is IPA:
% Ref: The Index of Pupillary Activity: Measuring Cognitive Load vis-à-vis Task Difficulty with Pupil Oscillation
% https://dl.acm.org/citation.cfm?id=3173856
% It has an open-source algorithm in Python which could be used instead of
% this