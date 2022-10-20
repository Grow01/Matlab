function [newdata] = lineardetrend(data,N,n_fin)

Nk=N/n_fin;
newmatrix = reshape(data,[Nk,n_fin]);   
D = detrend(newmatrix);  
newdata = reshape(D,[1,N]);

end