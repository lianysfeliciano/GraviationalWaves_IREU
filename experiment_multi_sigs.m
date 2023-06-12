function [p_times,p_freqs,p_eq_ps,xgrid,kgrid,color_grid,found] = experiment_multi_sigs...
    (N, minf,maxf,mc,white_noise,amp,plot_flag,save_fig,where_am_i )
%example 
% The code creates a source, simulates it, creates a peakmap, runs the
% Hough, selects candidates, and does coincidences with injection and
% candidates. This is used to create the plots in the paper with minf=140,
% maxf=150, mc=1e-3, amp=1e-22 (plots= peakmap, Hough map)

%   -----------------------------INPUTS:-----------------------------------
%   minf:           minimum frequency                                  (Hz)
%   maxf:           maximum frequency                                  (Hz)
%   mc:             chirp mass                                       (Msun)
%   white_noise     White noise (1)
%   amp:            amplitude of signal
%   plot_flag:      do plots (=1), else no
%   save_fig:       save sim figures
%   where_am_i:     cluster on which you run ('laptop','louvain','cnaf')

%   -----------------------------OUTPUTS:----------------------------------

%   n_sigs          number of retrieved signals
%   TTFT            
%   Plot 
%
%   -----------------------------EXAMPLE:----------------------------------

% example()

if ~exist('N','var')
    N = 1;
end

if ~exist('minf','var')
    minf=140;
end
if ~exist('maxf','var')
    maxf=150;
end
if ~exist('mc','var')
    mc=1e-3;
end
if ~exist('amp','var')
    amp=1e-22;
end

if ~exist('which_hough','var')
    which_hough='gfh';
end

if ~exist('plot_flag','var')
    plot_flag=1;
end

if ~exist('white_noise','var')
    white_noise = 0;
end

if ~exist('save_fig','var')                                    % option to save peak map and hough mapp
   save_fig=0;                         % Directory where figs are stored
end

if ~exist('where_am_i','var')
    where_am_i='laptop';
end
fap = 0.01;

n=11/3;                                                                     

if strcmp(where_am_i,'laptop') % my laptop
    
    sfdbdir='~/Nikhef_REU/sfdbs/H/' ; % directory where SFDBs are  

end


t00=1238777856+12000;                                                       % start time of injection
coin_dist=3;                                                                % distance in bins allowed between injection parms and cand parms
% rootkcand=10;                                                                % sqrt(num of candidates to select)
ref_perc_time=0;                                                          % reference time to use when constructing hough map: f0 is at t0


fdotmin=calc_fdot_chirp(mc,minf);                                           % fdot at minf and maxf
fdotmax=calc_fdot_chirp(mc,maxf);
TFFTmax=floor(1/sqrt(fdotmax));                                             % calculate TFFT using maximum spinup
TFFT=2^floor(log2(TFFTmax))                                                 % will work for TFFT not a factor of 2

t1=calc_time_to_coalescence(mc,minf);                                       % seconds until coalescence at min and max frequency injected
t2=calc_time_to_coalescence(mc,maxf);
dur=floor(t1-t2);                                                           % signal duration to simulate (s)


t0=gps2mjd(t00);                                                            % converts t00 (GPS time) to MJD time
sour=gen_N_power_law_sigs(N,minf,maxf,fdotmin,fdotmin,n,n,amp,t0);          % creates source structure, for chirp sigs, sig's f0 alwyas at minf

% for ss = 1:length(sour)
%     tt = 0:dur;
%     f_after_dur = power_law(tt,sour(ss).f0,sour(ss).df0,sour(ss).n);
%     if f_after_dur(end) >maxf
%         [~,b] = min(abs(f_after_dur-maxf));
%         tgw = tt(b);
%         ffdoot(ss) = calc_fdot_chirp(sour(ss).mc,f_max);
%     end
% end

% sour.dur=dur;

% sour(1).f0 = sour(1).f0+1/N;                                                  % to remove from the edge
% sour(1).x0 = 1/sour(1).f0^(n-1)

sig_type='power_law';

if strcmp(which_hough,'true_GFH')
%     thr = 0;
%     peak_flag = 0;
%     zero_flag = 1;
%     job_pack_0=inject_power_law_signal...
%         (sour,[minf maxf],t00,dur,TFFT,sfdbdir,sig_type,white_noise,thr,peak_flag,zero_flag);  
else
job_pack_0=inject_power_law_signal...
    (sour,[minf maxf],t00,dur,TFFT,sfdbdir,sig_type,white_noise);           % does the injection of sour in freq band, makes peakmap with TFFT;
end    
% job_pack_0 contains the peakmap
basic_info=job_pack_0.basic_info;                                           % structure that contains basic information about the peakmap    


% [SLong, SLat]=astro_coord('equ','ecl',sour.a,sour.d);                       % conversation of equatorial to ecliptic sky coordinates

p=job_pack_0.peaks;                                                         % the peaks in time/frequency, with corresponding amps
p_times=p(1,:);
p_freqs=p(2,:);
p_eq_ps=p(3,:);
if plot_flag==1
    plot_triplets(86400*(p(1,:)-p(1,1)),p(2,:),p(3,:))
    xlabel('time (s)'); ylabel('frequency (Hz)'); cblabel('equalized spectrum')
    set(gca,'FontSize',14)
end
pout=p;%andrew_hfdf_patch(p,basic_info,[SLong SLat]);                          % Doppler correction

hm_job.minf=minf;                                                           % minimum frequency to do the Hough on
hm_job.maxf=maxf;                                                           % maximum frequency to do the Hough on
hm_job.df=1/TFFT;                                                           % step in frequency
hm_job.dur=dur;
try
    hm_job.patch=[SLong SLat];
catch
    hm_job.patch = [0 0];
end
hm_job.n=n;
hm_job.ref_perc_time=ref_perc_time;

if strcmp(which_hough,'gfh')
    hm_job.frenh=1;
    [gridk,~]=andrew_long_transient_grid_k(TFFT,[minf maxf],[fdotmin fdotmax],dur,n);

    if n == 11/3
        kss_inj = [sour.kn]; mink = min(kss_inj); maxk = max(kss_inj);
        gridk = cbc_shorten_gridk(gridk,mink,maxk);
    end
    if length(gridk)>10000 
        gridk=shorten_gridk(gridk,sour.kn);                                 % restricts gridk to values around k of source
    end

    hm_job.gridk=gridk;
    tic;
    [hfdf,hm_job]=hfdf_hough_transients(pout,hm_job);                       % does the Generalized Hough with nonuniform k grid, uniform x grid
    toc;
    Nparmspace = numel(y_gd2(hfdf));
    rootkcand = ceil(sqrt(fap * Nparmspace));                                          % sqrt(num of candidates to select)                                                                        
    [cand2, job_info]=hfdf_peak_transients(hfdf,hm_job,rootkcand,0,0);      % selects candidates in every box of hough map of size rootkcand x rootkcand
    cand_orig=cand2;
  






    %%

    for ss = 1:length(sour)
    
        cand2=shift_pl_cand_parms(cand2,sour(ss),ref_perc_time);                    % shifts candidate x0 / f0 to reference time of inj, which is at t=0
        if sour(ss).n==11/3
            sour(ss).kn=-sour(ss).kn;                                                   % flip sign of k so that when doing coincidences with injections, things work
        end
        [found(ss),distaway(ss),dist_each_parm(ss,:),best_cand(:,ss)]=quick_coins(cand2,job_info,sour(ss),coin_dist); % coincidences between source and each candidate
        found(ss)
        
        x0(ss)=best_cand(1,ss);
        t0(ss)=best_cand(9,ss);
        kn(ss)=best_cand(4,ss);
        f0(ss)=best_cand(14,ss);
% best_cand_orig(:,ss)=cand_orig(:,x0(ss)==cand2(1,:));
%fdot0=best_cand(15,1);
        fdot0(ss)=abs(kn(ss)*f0(ss)^best_cand(7,ss));
        mc_found(ss)=calc_mc_with_k(abs(kn(ss)));
    end
    
    xgrid=x_gd2(hfdf);
    kgrid=gridk;
    color_grid=y_gd2(hfdf);
    
elseif strcmp(which_hough,'fh')
    hm_job.oper='noadapt';                                                  % The Hough is nonadaptive, meaning that the antenna pattern is not considered
    kcand=20;                                                               % 20 candidates are selected from the Frequency-Hough map uniformaly in f and fdot
    hm_job.frenh=10;                                                        % Frequency-Hough uses overesolution factor of frequency to improve sensitivity
    hm_job.sd=construct_grid_sd_parms(sour.df0);                            % Creates spin-up or spin-down grid for Frequency-Hough
    [hfdf,~,~,hm_job]=fu_transients_hough(pout,hm_job,ref_perc_time);       % Does the Frequency-Hough
    [cand2, job_info]=fu_transients_hfdf_peak(hfdf,hm_job, kcand);          % selects candidate in Frequency-Hough map

 
    %%% shifts the candidate to t=0 reference time because sour is at t=0
    [ cand2 ] = shift_cand_parms( cand2,sour,ref_perc_time );               % shifts candidate f0 to reference time of inj, which is at t=0
    [ found,best_cand,false_alarm,distaway,dist,dist_each_parm_best_inj ] = fu_coins_inj(cand2,job_info,sour,coin_dist);  % coincidences between source and each candidate
    found

elseif strcmp(which_hough,'true_GFH')
  
end




%%Creating Peakmap and Hough Map                                                                    

if plot_flag==1
    peak_map=plot_triplets(86400*(p(1,:)-p(1,1)),p(2,:),p(3,:))                                  %plots peakmap | frequency vs time
    xlabel('time (s)'); ylabel('frequency (Hz)'); cblabel('equalized spectrum')
    set(gca,'FontSize',14)
end

if plot_flag==1
   if strcmp(which_hough,'gfh')
        hough_map = image_gd2_nonunigrid(hfdf,calc_mc_with_k(gridk),n);                           % plots the hough map, k vs. f, number count colored
        xlabel('frequency (Hz)'); ylabel('chirp mass (M_{sun})'); cblabel('number count')
   end
end



%%Saveing Figures
if save_fig==1
    saveas(peak_map,'/Output_Code/Peak_map.png')
    saveas(hough_map,'/Output_Code/Peak_map.png')
end



end

