function [df] =Sim_df(Num_sims, n,Retireved_sig,tfft,minfreq,maxfreq,amps,chirpm)

%Sim_df
%This function creates and stores a dafa frame of a hough simulation
%   This function is intended to store values of many simulations at once. 
%   By implimenting this function within a loop of simulations it will
%   store every loop
% 
%-----------------------------------------------------------------------
%               
% df                   Dataframe with all input values
%
%df is separated by commas and has a header of one row. No idx col. 
%-------------------------------------------------------------------------------


%Creating empty arrays for each column

 N=zeros(Num_sims,1);
 Retireved_N=zeros(Num_sims,1); 
 TFFT=zeros(Num_sims,1);
 minf=zeros(Num_sims,1);
 maxf=zeros(Num_sims,1);
 amp=zeros(Num_sims,1);
 cm=zeros(Num_sims,1);




%Insertinng values into empty arrays
for i = 1:Num_sims
    N(i)=n(i);
    Retireved_N(i)=Retireved_sig(i); 
    TFFT(i)=tfft(i);
    minf(i)= minfreq(i);
    maxf(i)=maxfreq(i); 
    amp(i)=amps(i);
    cm(i)=chirpm(i);


end


%Composing tabel of values
df=table(N,Retireved_N,TFFT,minf,maxf,amp,cm);

%Saving df as csv file 
writetable(df,'df.csv');

end
