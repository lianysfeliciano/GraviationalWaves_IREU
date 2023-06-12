%%
path=genpath('/Users/lianysfeliciano/Nikhef_REU');
addpath(path)

%[p_times,p_freqs,p_eq_ps,xgrid,kgrid,color_grid,found] = experiment_multi_sigs...
    %(N, minf,maxf,mc,white_noise,amp,plot_flag,save_fig,where_am_i )

 
%
i=1; 
p_times =zeros(4);
p_freqs =zeros(4);
p_eq_ps =zeros(4);
xgrind =zeros(4);
color_grind =zeros(4);
found =zeros(4);

while i<4  
    
    [times,freqs,p_eq_ps,xgrid,kgrid,color_grid,found]=example_multi_sigs(2,6,20,1.5,1,1e-22,1,1,'laptop');
    i=i+1;
    p_timee{i}=times;
    p_freqs{i}=freqs;
    xgrind{i}=xgrind;
    kgrind{i}=kgrind;
    color_grind{i}=color_grind
    found{i}=found

end 
