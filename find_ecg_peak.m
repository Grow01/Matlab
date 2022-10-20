function [values,Maxindex] = find_ecg_peak(newdata,S)

[values,Maxindex]=findpeaks(newdata,'MinPeakHeight',S);

end