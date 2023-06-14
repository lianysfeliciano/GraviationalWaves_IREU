%%
path=genpath('/Users/lianysfeliciano/Nikhef_REU');
addpath(path)

%[p_times,p_freqs,p_eq_ps,xgrid,kgrid,color_grid,found] = experiment_multi_sigs...
    %(N, minf,maxf,mc,white_noise,amp,plot_flag,save_fig,where_am_i )

%for i=1:4  
    
    %example_multi_sigs(); 

%end 

%Input Parameters
Num =1:10;
whitenoise=1;
minfreq=zeros(10,1); %Making these empty arrays to fill with all the same value | When varyinng these will be input arrays 
maxfreq=zeros(10,1);
chirpm=zeros(10,1);
amplitute=zeros(10,1);
donde='laptop';

%Filling Empty Input Arrays

minfreq(:)=3;
maxfreq(:)=19;
chirpm(:)=2.5;
amplitute(:)=1e-25;


%Creating empty arrays for outputs
registered_sig= zeros(10,1);
%tfft= zeros(10,1);


for i=1:length(Num)
    
    %Simulating
    [p_times,p_freqs,p_eq_ps,xgrid,kgrid,color_grid,found]=experiment_multi_sigs(Num(i),minfreq,maxfreq,chirpm,whitenoise,amplitute,1,0,donde);
    
    %Filling output array with output vals
    registered_sig(i)=found;
    % tfft(i)=TFFTmax;
    
end

%Storing everything into a df
Sim_df(10,Num,registered_sig,tfft,minffreq,maxfreq,amplitute,chirpm)


