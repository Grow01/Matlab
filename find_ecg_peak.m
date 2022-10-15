function [pks, loc_pks]= find_ecg_peak(data)
    soglia=find_thr(data,6);
    [pks,loc_pks] = findpeaks(data,"MinPeakHeight",soglia);
end

function thr=find_thr(data, multiplier)
    media=mean(data);
    thr=media+multiplier*std(data);
end
