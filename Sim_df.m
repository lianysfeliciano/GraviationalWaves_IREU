function [df] =Sim_df(Num_sims,n,Retireved_sig,tfft,Fo,Dist,Durration,cr,chirpm,minfreq,maxfreq,amps)

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
 TFFT=zeros(Num_sims,1);               %lengths of time of Foyer transmform 
 minf=zeros(Num_sims,1);               %sim input
 maxf=zeros(Num_sims,1);               %sim input
 h0=zeros(Num_sims,1);                 %amplitude
 cm=zeros(Num_sims,1);                 %sim input
 f0=zeros(Num_sims,1);                 %frequency at given points (this will be a nested array in each entry)
 dur=zeros(Num_sims,1);                % durration of time signal is ditected
 Cr=zeros(Num_sims,1);                 %critical ratio
 %distaway=zeros(Num_sims,1);

%Insertinng values into empty arrays
for i = 1:Num_sims
    N(:)=n;
    Retireved_N(i)=Retireved_sig(i); 
    TFFT(i)=tfft(i);
    minf(i)= minfreq(i);
    maxf(i)=maxfreq(i);
    Cr(i)=cr(i);
    cm(i)=chirpm(i);
    h0(i)=amps(i);
    dur(i)=Durration(i);
    %f0(i)=Fo(i);  
    %distaway(i)=Dist(i);

end


%Composing tabel of values

df=table(N,Retireved_N,TFFT,dur,Cr,cm,minf,maxf,h0);

%Merging tables into one large dataframe
df=[df,Dist];
df=[df,Fo];


N1=N(1);
h01=h0(1);

%Naming the datafile
filename= append('df_',num2str(N1),"Sig_Amp:",num2str(h01));

%Saving df as csv file 
writetable(df,append(filename,'.csv'));


end






