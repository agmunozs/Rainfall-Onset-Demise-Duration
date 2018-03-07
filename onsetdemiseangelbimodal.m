% Matlab script to compute onset/demise (and duration)
% AG Munoz (IRI Columbia U) - agmunoz@iri.columbia.edu and Carlos Martinez - carlos.martinez@columbia.edu
% Project: Predictability of the onset, duration and demise of the North and South American Monsoon Systems: 
% the role of cross-equatorial interactions
% First edition: Aug 4, 2017 
% Last edition: Feb 21, 2018
% Notes: 
% + Subset of the full code written for the N/S American Monsoon
%   project, to show the Bombardi and Carvalho approach. 
% + Code reads rainfall data for indicated coordinates (1 gridbox) directly
%   from the IRI Data Library (CPC-Unified - Chen et al 2008).
% + For illustration purposes, only one gridbox/station is used.
% + The method has issues when the year is too dry (as one should expect).
%3,4,14,16,17,18,19,20
%%%%%START OF USER-MODIFIABLE SECTION%%%%%%%%%%%%
% 38 Stations in the Caribbean
% 57 Years of Data (1960 -  2017)
statn = 25;% Which station I want
sf=10; %smooth factor for the sgolay smoother
names = {'Antigua_Barbuda_IntlAP' ,'Nassau_IntlAP_Bahamas' ,'CIMH_Barbados', 'Gantley_Barbados', 'C_Farm_Belize',  'Intl_AP_Belize', 'Georgetown_Caymen'...
    , 'Camaguey_Cuba', 'La_Habana_Cuba', 'DCAP_Dominica', 'Santo_Domingo_DR', 'Guadeloupe_IntlAP_Guadeloupe', 'Georgetown_Guyana', 'Time_HRI_Guyana' ...
    , 'Worthy_Park_Jamaica', 'IntlAP_Martinique', 'Henaworra_IntlAP_St_Lucia', 'Dumbarton_St_Vincent', 'Corantjinpolder_Suriname', 'Zanderji_Suriname', 'Piarco_IntlAP_TT' ...
    , 'Crown_Point_TT', 'Ft_Lauderdale_USA', 'Key_West_USA', 'Miami_IntlAP_USA', 'Palm_Beach_USA', 'Henry_E_R_IntilAP_St_Criox', 'Cyril_E_King_IntlAP_St_Thomas' ...
    , 'Colosso_USPR', 'Dora_Bora_USPR', 'Ensenda_USPR', 'Guaynama_USPR', 'JaJome_Alto_USPR','Mora_Camp_USPR','Paraiso_USPR' ...
    , 'Morovis_USPR', 'San_Andreas_Colombia', 'Felipe_Mexico'};

load dailydata.mat
%%%END OF USER-MODIFIABLE SECTION (DO NOT MODIFY ANYTHING BELOW THIS LINE)%%%%%
%%
%Daily and Pentad Record
%Reshape for Daily Data into 365 days for 57 years
for j = 1:38
yryrdaily.(names{j}) = reshape(dailydata(:,j), 365, []); 
end

%Pentad Means
for i = 1:73
    for k = 1:57
        for j = 1:38
           if nansum(isnan(yryrdaily.(names{j})(1+5*(i-1):5+5*(i-1),k))) > 2 %If pentad has 3 or more NaN's, then pentad is NaN
            yrypentadmean.(names{j})(i,k) = NaN;
           else
            yrypentadmean.(names{j})(i,k) = nanmean(yryrdaily.(names{j})(1+5*(i-1):5+5*(i-1),k));
           end
        end
    end
end

%Pentad Climatology
for i = 1:73
    for j = 1:38
        pentclim(i,j) = nanmean(yrypentadmean.(names{j})(i,:));
    end
end
% pentclim replaces what was used as penmean
% pentad is replaced by the yrypentadmean

%% EDITED 3/6/2018: Annual Climatological Precipitation Rate for PC
for j = 1:38
    PC(j) = mean(pentclim(:,j));
end
%% Angel M. Onset and Demise Code
% %Computing onset/demise a la Bombardi and Carvalho 2009:
clear onset demise inipent inipent2
pentad = yrypentadmean.(names{statn})(:,:); % Calls in the 29th Station (Coloso, PR) 73x57yrs
ns=size(pentad);nyear=ns(2);

    
%Find initial pentad for onset
display('Computing initial pentad... (this takes some time)')

%inipent/inipent2 is t0 in Bombardi and Carvalho 2009, initial pentad for
%onset and demise, respectively. Here wrt Jan 1 every year

        clear f1 f2
            if(~(isnan(pentclim(:,statn))))
                %Compute first harmonics of the climatological rainfall for
                %the location:
                f1 = fit([1:73]',squeeze(pentclim(1:73,statn)),'fourier3');
                yfitted = feval(f1,[1:73]);
                [ypkon,idxon] = findpeaks(-yfitted);%min, for onset
                [ypkde,idxde] = findpeaks(yfitted);%max, for demise
                %Find x-coordinates of two larger maxima
                [so,al]=sort(ypkde); %two maxima are at the end of so; al is location
                %idxde(al(end))
                %idxde(al(end-1))
                %Find minimum before first large maximum, for rainfall onset
                xfirstmax=min(idxde(al(end)),idxde(al(end-1)));
                pp=find(idxon < xfirstmax);
                %% EDITED 3/6/2018 CHANGED NAME TO INIPENT0/02
                %Created two t0's for each inipent: One that begins when t0
                % = 2 in order to see the entire S,Sb,dSb time series, and
                % t0 = minima or maxima known as inipent0/02 that is used
                % ONLY for the ONSET/DEMISE Section.
                
                
                inipent0=idxon(pp(end));
                %Find maximum that is second larg maximum, for rainfall demise
                xsecmax=max(idxde(al(end)),idxde(al(end-1)));
                %pp=find(idxon > xsecmax);
                inipent02 = xsecmax;
                %inipent2=idxon(pp(1));
                %Find minimum between two largest maxima, for MSD purposes
                pp=find(idxon < xsecmax & idxon > xfirstmax);
                inipent3=idxon(pp);
     
                %----
                %Older versions:
                %inipent=round(fminbnd(f1,1,73)); %min, for onset
                %inipent=round(interp1(feval(f1,[1:73]'),1:73,min(feval(f1,1:73))));%min, for onset
                %f2 = fit([1:73]',squeeze(-penmean(1:73)),'fourier2');  %now negative for maximum
                %inipent2=round(fminbnd(f2,1,73));%max, for demise
                %inipent2=round(interp1(feval(f1,[1:73]'),1:73,max(feval(f1,1:73))));%max, for demise
                %----
            else
                inipent=NaN;
                inipent2=NaN;
                inipent3=NaN;
            end
            
inipent0(inipent0==0)=NaN;
inipent02(inipent02==0)=NaN;
display('t0 Done!')
figure(1);clf %Changed mean to nanmean for (pentad,2)
plot(f1);hold on; plot(1:73,nanmean(pentad,2))%,inipent,f1(inipent),'*k',inipent2,f1(inipent2),'*k','markers',12)

                findpeaks(yfitted)
                findpeaks(-yfitted)

%%
%Compute S
clear S S0
display('Computing S...')
%% 3/6/2018 CHANGE WHERE INIPENT = 1 IN ORDER TO SEE ENTIRE S,SB TIMESERIES.
inipent = 1;
inipent2= 1;
        for ity=1:nyear-1 
            S0=0;
            S02=0;
            %if(isnan(penmean) | isnan(inipent) | isnan(inipent2))
            %    S(ity)=NaN;
            %    S2(ity)=NaN;   %S2 for demise, using inipent2 
            %else
                for itp=1:73  %counter for 
                    itk=inipent+itp;
                    itk2=inipent2+itp;  %for demise, using inipent2 
                    if itp == 1
                        %% 3/6/2018 CHANGE -- PC NOT PENTCLIM
                        S(ity,inipent)=squeeze(pentad(inipent,ity)-PC(statn)); 
                        S2(ity,inipent2)=pentad(inipent2,ity)-PC(statn);
                    else
                        if itk<=73
                            S(ity,inipent+itp-1)=nansum(pentad(inipent:itk,ity)-PC(statn),1);
                        else %use next year's data  
                            S0=nansum(pentad(inipent:73,ity)-PC(statn),1); %accumulated until the end of previous year
                            S(ity,inipent+itp-1)= S0 + (nansum(pentad(1:itk-73,ity+1)-PC(statn),1) ); %total, considering demise on following year
                         end
                        if itk2<=73  %demise
                            S2(ity,inipent2+itp-1)=nansum(pentad(inipent2:itk2,ity)-PC(statn),1);
                        else %use next year's data  
                            S02=nansum(pentad(inipent2:73,ity)-PC(statn),1); %accumulated until the end of previous year
                            S2(ity,inipent2+itp-1)= S02 + (nansum(pentad(1:itk2-73,ity+1)-PC(statn),1) ); %total, considering demise on following year
                        end
                    end
                end
            %end
        end

S(S==0)=NaN;
S2(S2==0)=NaN;
display('S Done!')

%Compute Sb
clear Sb Sb2
display('Smoothing with a Savitzky-Golay 2nd-degree filter...')
%% 3/6/2018 CHANGE WHERE I MAKE INIPENT RANGE IN S FROM 1:73 -- CAN BE ALTERED TO FOCUS ON CERTAIN RANGES.
        for ity=1:nyear-1
            if( isnan(S(ity,inipent:72+inipent)) | isnan(S2(ity,inipent2:72+inipent2)) )
                Sb(ity,:)=NaN;
                Sb2(ity,:)=NaN;
            else
                %for pentads, R. Bombardi suggests to smooth S ~20 times;
                %for daily data ~150 times !!! NOT USING THAT, USING SGOLAY
                 Sb(ity,inipent:72+inipent)=smooth(S(ity,inipent:72+inipent),sf,'sgolay');
                 Sb2(ity,inipent2:72+inipent2)=smooth(S2(ity,inipent2:72+inipent2),sf,'sgolay');
            end
        end

display('Sb Done!')

%Compute dSb/dt
clear dSb 
%zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);    %function to find zero-crossings   --NOT USED ANYMORE

        for ity=1:nyear-1 
            if(isnan(Sb(ity,inipent:72+inipent)))
                dSb(ity,:)=NaN;
            else
                dSb(ity,inipent:71+inipent)=diff(Sb(ity,inipent:72+inipent));   
            end
            if(isnan(Sb2(ity,inipent2:72+inipent2)))
                dSb2(ity,:)=NaN;
            else
                dSb2(ity,inipent2:71+inipent2)=diff(Sb2(ity,inipent2:72+inipent2));   
            end
        end

display('dSb Done!')



%% 3/6/2018 CHANGE: Investigates dSb after starting dates inipent0 and 02
% This allows to see where the onset and demise is AFTER t0. (filter out)
dSbomit = zeros(size(dSb));
dSbomit(:,inipent0-1:end) = dSb(:,inipent0-1:end);

dSbomit2 = zeros(size(dSb2));
dSbomit2(:,inipent02-1:end) = dSb2(:,inipent02-1:end);

%%%Compute onset and demise
clear onset demise onset4anom demise4anom
        for ity=1:nyear-1 
            o=0;
            d=0;
            if(isnan(dSb(ity,inipent:71)) | isnan(dSb2(ity,inipent2:71)))
                onset(ity)=NaN;
                demise(ity)=NaN;
            else
                for itp = inipent0:71 %itp=2:71  %pentad relative to t0!!!
                    if o==1
                    else
                        if (dSbomit(ity,itp-1) < 0 && dSbomit(ity,itp) > 0 && dSbomit(ity,itp+1) > 0 && dSbomit(ity,itp+2) > 0 )
                        onset(ity)=itp+1;   %% RELATIVE TO t0
                        onset4anom=onset;
                          if onset(ity)>73
                              onset(ity)=onset(ity)-73;
                          end
                        o=1;
                        else
                          onset(ity)=NaN;  
                        end
                    end
                end
            end
                 for itp = inipent02:71 %% ADDED THIS IN ORDER TO CHANGE RANGES IF ONSET/DEMISE IS NOT REASONABLY RIGHT.
                    if d==1
                    else
                 
                        if ( dSbomit2(ity,itp-1) > 0 && dSbomit2(ity,itp) < 0 && dSbomit2(ity,itp+1) < 0  && dSbomit2(ity,itp+2) < 0)
                        demise(ity)=itp-1;  %relative to Jan 1st!
                        demise4anom=demise;
                        %if ( dSb(ilo,ila,ity,itp-1) > 0 && dSb(ilo,ila,ity,itp) < 0  && dSb(ilo,ila,ity,itp+1) < 0 )
                        %demise(ilo,ila,ity)=itp+inipent(ilo,ila);  %relative to Jan 1st!
                          if demise(ity)>73
                              demise(ity)=demise(ity)-73;
                          end
                        d=1;
                        else
                          demise(ity)=NaN;
                        end
                    end
                end
        end

display('Onset and Demise Done!')
save -v7.3 onsetdemise.mat onset onset4anom demise demise4anom inipent inipent2 S Sb dSb pentad 

%%
%TEST PLOTS
%Select a year:
iyr=53;
%
figure(7);clf %S,Sb,dSb for onset
fig = figure(7);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
bfig=bar(1:72,squeeze(pentad(1:72,iyr)')*1, .8, 'Facecolor',colors('celestial blue'), 'Edgecolor',colors('black')); hold on
plot(onset(iyr),1,'*g',demise(iyr),1,'*k','markers',12,'LineWidth', 2), hold on;
ylabel('Precipitation (mm/pentad)')
yyaxis right
plot(1:72,squeeze(S(iyr,1:72)),'k', 'LineWidth', 3), hold on;
plot(1:72,squeeze(Sb(iyr,1:72)),'Color',[.91,.41,.17],'LineWidth',3);
plot(1:72,squeeze(dSb(iyr,1:72)),'Color',[0.05,.7,.4],'LineWidth',3, 'Linestyle','-'), hold on;
line([1,73],[0,0], 'color', colors('burgundy'), 'LineStyle', '--', 'LineWidth', 3);
line([onset(iyr),onset(iyr)],[-100,0], 'color', colors('burgundy'), 'LineStyle', '--', 'LineWidth', 3);
line([demise(iyr),demise(iyr)],[-100,0], 'color', colors('burgundy'), 'LineStyle', '--', 'LineWidth', 3);
ylabel('S,Sb,dSb for onset')
title(names(statn))
xlabel('pentad')
legend('rainfall','onset','demise','S','Sb','dSb/dt','Zero Line','Onset Line','Demise Line')
legend('Location','northwest')
legend('boxoff')

% figure(3);clf %S,Sb,dSb for demise
% bfig=bar(1:72,squeeze(pentad(1:72,iyr)')*1);   
% ylabel('Precipitation (mm/pentad)')
% yyaxis right
% pfig=plot(1:72,squeeze(S2(iyr,1:72)),1:72,squeeze(Sb2(iyr,1:72)),1:72,squeeze(dSb2(iyr,1:72)),'LineWidth',3);hold on
% ylabel('S,Sb,dSb for demise') 
% plot(onset(iyr),0,'*g',demise(iyr),0,'*k','markers',15,'LineWidth', 2)
% title(names(statn))
%  xlabel('pentad')
%  legend('S','Sb','dSb/dt','zero line','onset','demise')
%  legend('Location','northeast')
%  legend('boxoff')