function [df] =Sim_df(Name,Loc, N,Retireved_sig,Params,TFFT)

%Sim_df
%This function creates and stores a dafa frame of a hough simulation
%   This function is intended to store values of many simulations at once. 
%   By implimenting this function within a loop of simulations it will
%   store every indivisual loop
%  
%Inputs
%------------------------------------------------------------------------
% Variable Name        Definition
% 
% Name                 Name of Dataframe
% Loc                  Desiered Location of Dataframe
% Num_sims             Number of simulations run 
%  N                   Number of injected signals
% Retrived_sig         Number of retireved signals 
% Params               Array of Input Parameters (minf,maxf,Mc,white_noise, )
% TFFT                 Total number of TFFT's
% 
%
% -----------------------------------------------------------------------------
%
%Outputs
%-----------------------------------------------------------------------
% Variable Name        Definition
% df                   Dataframe with all input values
%
%-------------------------------------------------------------------------------

%Setting default values of input if none are provided. 

Name= ("/df_"+ N +"_signals");

Loc= "~/Nikhef/Output_Code";


%Creating empty arrays for each column
if N_sims>1
    N=zeros(Num_sims);
    Retireved_sig=zeros(Num_sims);
    TFFT=zeros(Num_sims);
end


%Insertinng values into 
for i = 1:N_sims
    N{i}=
    Retireved_sig{i}=
    TFFT{i}=

end

df=tabel(N,Retireved_sig,Params,TFFT);

saveas(df,Loc+Name+".csvread")






end