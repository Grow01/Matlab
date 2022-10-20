function [S] = find_thr(newdata,multiplier)

S=mean(newdata)+multiplier*std(newdata);
end