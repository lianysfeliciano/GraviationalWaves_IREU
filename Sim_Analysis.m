%%
path=genpath('/Users/lianysfeliciano/Nikhef_REU');
addpath(path)

%[p_times,p_freqs,p_eq_ps,xgrid,kgrid,color_grid,found] = experiment_multi_sigs...
    %(N, minf,maxf,mc,white_noise,amp,plot_flag,save_fig,where_am_i )
    
%for i=1:4  
    
    %example_multi_sigs(); 

%end 

%Number of Sims
Sims=50; 

%Innput Parameters
 
%For random values withinn range use this: xmin+rand(1,n)*(xmax-xmin)
Num =[1,2,3,4,5,6,7,8,9,10];
whitenoise=1;
minfreq=zeros(Sims,1);
maxfreq=zeros(Sims,1);
chirpm=zeros(Sims,1);
amplitute=[1e-22,1e-23,1e-24];


%minfreq=randi([2,4],Sims,1);      % Array of random (integer) values between 2 and 4 for min frequency (Must be whole numbers
%maxfreq=randi([6,8],Sims,1);      
%chirpm=0.5+rand(Sims,1)*(3-0.5);

donde='laptop';

%Filling Empty Input Arrays
minfreq(:)=4;
maxfreq(:)=7;
chirpm(:)= 1.15;


%Creating empty arrays for outputs
registered_sig= zeros(Sims,1);
tfft= zeros(Sims,1);
Dur=zeros(Sims,1);
Cr=zeros(Sims,1);
DistAway=zeros(Sims,1);
F0=zeros(Sims,1);
amp=zeros(Sims,1);

%Initialiting tables 



for n=1:10  
    
    for a=1:3
    Dist=table();
    F0=table();

        for i=1:Sims
        %Simulating
            [TFFT,dur,best_cand,distaway,found,f0]=experiment_multi_sigs(Num(n),minfreq(i),maxfreq(i),chirpm(i),whitenoise,amplitute(a),0,0,donde);
    
        %Filling output array with output vals
            registered_sig(i)=sum(found);
            tfft(i)=TFFT;
            Dur(i)=dur;
            Cr(i)=best_cand(6);
            %DistAway(i)=distaway;
            %F0(i)=f0;
            amp(i)=amplitute(a);
           
            
            %Table Variables

            f0_temp=table(f0);

            F0=[F0;f0_temp];

            dist_temp = table(distaway) ;
            Dist = [Dist ; dist_temp];
   
        end
        %Storing everything into a df
        Sim_df(Sims,Num(n),registered_sig,tfft,F0,Dist,Dur,Cr,chirpm,minfreq,maxfreq,amp)

    end

end






