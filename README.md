# Ocular_data_analysis

Ocular event detection for monocular video-based eye-trackers.

Here you can find MATLAB codes to detect saccades, fixations, blinks, and processing of the time-series of pupil diameter.

The codes are based on well-known algorithms in the literature, which are mentioned in-line with the codes.

You may need to customize the codes based on the eye-tracker you are using.

# Usage

You need to run ET_processing.m which loads the .m data called ET_data for example which is structured as follows:

The ocular data is assuemed to be converted into .mat file, e.g. with the name: ET_data 

It is also assumed to containt the following data arrays (vectors):

EH_gaze_hcoord: Horizontal Gaze Coordinate over time

EH_gaze_vcoord: Vertical Gaze Coordinate over time

EH_gaze_length: Perpendicular distance between the calibration surface and the participant/user over time

PplDiam: Pupil diameter over time (zero and NaN mean no detection of the pupil or closed eyes)

cr_diam/cr_det: Corneal diameter over time (zero and NaN mean no detection of the cornea or closed eyes)

scen_num: scene number over time (optional: to detect whether the gaze data is refering to the calibrated surface or not. It is zero for the default scene defined in ASL Eye-trac7)

The output variables include fixation samples (f_samples), blink samples (b_samples), and saccade samples (s_smpl). The calculation for some oculometrics are also available, e.g. saccade peak velocity (s_pv).

You need to specify the sampling frequency (Fs) in the script. (Default is 360 Hz)

# Further considerations:

Note that some of the criteria and thresholds may be modified depending on the task and settings of the recordings.

The interpolations can be used for visualization but not for the computation of oculometrics.

This algorithm parse the gaze data into saccades, fixations, blinks, and undefined events. 

It also provides pre-processings for pupil size.

Post-saccade oscillations are merged into fixations if they satisfy fixation criteria.

The detection of the ocular events is a sequential, meaning that the fixation detection is dependent to saccade detection.

The ocular data is assumed to be sampled regularly and corrected for possible drifts.

# References

Marandi, R.Z., Madeleine, P., Omland, Ø., Vuillerme, N. and Samani, A., 2018. Reliability of oculometrics during a mentally demanding task in young and old adults. IEEE Access, 6, pp.17500-17517.

Marandi, R.Z., Madeleine, P., Omland, Ø., Vuillerme, N. and Samani, A., 2019. An oculometrics-based biofeedback system to impede fatigue development during computer work: A proof-of-concept study. PloS one, 14(5).

Marandi, R.Z., Madeleine, P., Omland, Ø., Vuillerme, N. and Samani, A., 2018. Eye movement characteristics reflected fatigue development in both young and elderly individuals. Scientific reports, 8(1), p.13148.
