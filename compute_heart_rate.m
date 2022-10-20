function [bpm] = compute_heart_rate(Maxindex,tmax)

fs=length(Maxindex)/tmax;
bpm=fs*60;

end