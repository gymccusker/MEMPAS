%%
%
%       SIZE-SEGREGATION FUNCTION
%

function [OUTPUT,zz]=sem_sizeseg(OUTPUT,flag)

% % % % % % % % % % % % % %
% % % % % % % % % % % % % %
% % % % % % % % % % % % % %       Size-Segregated Composition Distributions
% % % % % % % % % % % % % %   
% % % % % % % % % % % % % %

% PCASP and CDP combined, PCASP small bins averaged
bins=[0 0.1302 0.1551 0.2480 0.4682 0.6347 0.8097 0.9155 1.0214 1.2474 1.4514 1.7141 1.8594 2.0087...
    2.1930 2.3431 2.5241 3.7100 4.6200 6.0000 7.3700 8.4900 9.0500 9.6500 10.7200 10.9200 10.9700 12.6500...
    15.5800 18.3500 19.7900 20.5200 22.1500 24.7400 27.3300 29.3500 31.0500 32.7700 34.6600 37.2200 39.5100...
    41.1600 42.8700 44.8400 46.6100 48.6600]; 

%   0.5um - 10um 
bins_new=[0 bins(6:25)];


bins_switch = input('All bins (1), or 0.5um to 10um only (2): ');
switch bins_switch
    case 1
        
        % % % % % %
        % % % % % %     CENTRE ON PROBE BINS
        % % % % % %
        
        case_label=1;
%         lim_binwidth=bins(2:end)-bins(1:end-1);
        lim_upper=[bins(1:end-1)+(bins(2:end)-bins(1:end-1))/2 bins(end)+(bins(end)-bins(end-1))/2];
        lim_lower=[bins(1) bins(2:end)-(bins(2:end)-bins(1:end-1))/2];
        
        u=cell(1,length(bins));
        for i=1:length(bins),
            u{1,i}=num2str(bins(i),'%.2f');
        end
        ll=u(1:4:end);
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                
    case 2

        % % % % % %
        % % % % % %     CENTRE ON PROBE BINS
        % % % % % %        
        case_label=2;
%         lim_binwidth=bins_new(2:end)-bins_new(1:end-1);
        lim_upper=[bins_new(1:end-1)+(bins_new(2:end)-bins_new(1:end-1))/2 bins_new(end)+(bins_new(end)-bins_new(end-1))/2];
        lim_lower=[bins_new(1) bins_new(2:end)-(bins_new(2:end)-bins_new(1:end-1))/2];
             
        w=cell(1,length(bins_new));
        for i=1:length(bins_new),
            w{1,i}=num2str(bins_new(i),'%.2f');
        end
        ll=w(1:2:end);
%         ll={'0.5' '0.8' '1.0' '1.5' '1.9' '2.2' '2.5' '4.6' '7.4' '9.0' '10.7'};

        
    otherwise
            disp('CHOICE NOT AVAILABLE')
end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        
        if flag==1,

        for i=1:length(lim_upper(1:end-1)),
            eval(['dust.ind',num2str(i,'%d'),'=find(OUTPUT.Silicate.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Silicate.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['ca.ind',num2str(i,'%d'),'=find(OUTPUT.CaRich.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.CaRich.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['freshcl.ind',num2str(i,'%d'),'=find(OUTPUT.FreshCl.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.FreshCl.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['agedcl.ind',num2str(i,'%d'),'=find(OUTPUT.AgedCl.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.AgedCl.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['sulph.ind',num2str(i,'%d'),'=find(OUTPUT.Sulph.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Sulph.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['mixedsil.ind',num2str(i,'%d'),'=find(OUTPUT.MixedSil.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.MixedSil.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['carbon.ind',num2str(i,'%d'),'=find(OUTPUT.Carbon.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Carbon.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['bio.ind',num2str(i,'%d'),'=find(OUTPUT.Bio.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Bio.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['phos.ind',num2str(i,'%d'),'=find(OUTPUT.Phos.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Phos.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['metal.ind',num2str(i,'%d'),'=find(OUTPUT.Metal.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Metal.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['gypsuml.ind',num2str(i,'%d'),'=find(OUTPUT.Gypsum.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Gypsum.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['biomass.ind',num2str(i,'%d'),'=find(OUTPUT.Biomass.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Biomass.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['other.ind',num2str(i,'%d'),'=find(OUTPUT.Other.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Other.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['totl.ind',num2str(i,'%d'),'=find(OUTPUT.All.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.All.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['ns',num2str(i,'%d'),'_bins=length(totl.ind',num2str(i,'%d'),');'])
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['inpdust_bins(',num2str(i,'%d'),')=length(dust.ind',num2str(i,'%d'),');']);
            eval(['inpmix_bins(',num2str(i,'%d'),')=length(mixedsil.ind',num2str(i,'%d'),');']);
            eval(['inpca_bins(',num2str(i,'%d'),')=length(ca.ind',num2str(i,'%d'),');']);
            eval(['inpfreshss_bins(',num2str(i,'%d'),')=length(freshcl.ind',num2str(i,'%d'),');']);
            eval(['inpagedss_bins(',num2str(i,'%d'),')=length(agedcl.ind',num2str(i,'%d'),');']);
            eval(['inpsulphate_bins(',num2str(i,'%d'),')=length(sulph.ind',num2str(i,'%d'),');']);
            eval(['inpgypsum_bins(',num2str(i,'%d'),')=length(gypsuml.ind',num2str(i,'%d'),');']);
            eval(['inpcarbon_bins(',num2str(i,'%d'),')=length(carbon.ind',num2str(i,'%d'),');']);
            eval(['inpbio_bins(',num2str(i,'%d'),')=length(bio.ind',num2str(i,'%d'),');']);
            eval(['inpphos_bins(',num2str(i,'%d'),')=length(phos.ind',num2str(i,'%d'),');']);
            eval(['inpox_bins(',num2str(i,'%d'),')=length(metal.ind',num2str(i,'%d'),');']);
            eval(['inpbiomass_bins(',num2str(i,'%d'),')=length(biomass.ind',num2str(i,'%d'),');']);
            eval(['inpother_bins(',num2str(i,'%d'),')=length(other.ind',num2str(i,'%d'),');']);
        end
        
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['inp',num2str(i,'%d'),'_bins=[inpdust_bins(',num2str(i,'%d'),') inpmix_bins(',num2str(i,'%d'),') inpca_bins(',num2str(i,'%d'),...
                ') inpfreshss_bins(',num2str(i,'%d'),') inpagedss_bins(',num2str(i,'%d'),') inpsulphate_bins(',num2str(i,'%d'),...
                ') inpgypsum_bins(',num2str(i,'%d'),') inpcarbon_bins(',num2str(i,'%d'),') inpbio_bins(',num2str(i,'%d'),...
                ') inpphos_bins(',num2str(i,'%d'),') inpox_bins(',num2str(i,'%d'),') inpbiomass_bins(',num2str(i,'%d'),...
                ') inpother_bins(',num2str(i,'%d'),')];']);
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['input',num2str(i,'%d'),'_bins=inp',num2str(i,'%d'),'_bins./sum(inp',...
                num2str(i,'%d'),'_bins);']);
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['mindustindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(1);']);
            eval(['mixindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(2);']);
            eval(['carichindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(3);']);
            eval(['freshssindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(4);']);
            eval(['agedssindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(5);']);
            eval(['sulphateindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(6);']);
            eval(['gypsumindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(7);']);
            eval(['carbonindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(8);']);
            eval(['bioindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(9);']);
            eval(['phosindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(10);']);
            eval(['oxindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(11);']);
            eval(['biomassindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(12);']);
            eval(['otherindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(13);']);
        end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  


zz=figure();
bar1=bar(1:length(lim_lower(1:end-1)),[mindustindex' mixindex' carichindex' freshssindex' agedssindex' sulphateindex' gypsumindex' ...
    carbonindex' bioindex' phosindex' oxindex' biomassindex' otherindex'],1.0,'stack');
if case_label==1, set(gca,'xtick',1:4:length(lim_lower)); end
if case_label==2; set(gca,'xtick',1:2:length(lim_lower)); end
set(gca,'XTickLabel',ll,'fontsize',16) 
set(bar1(1),'FaceColor',[0 0 1],'DisplayName','Silicates','EdgeColor',[0 0 0]);
set(bar1(2),...
    'FaceColor',[0.34901961684227 0.200000002980232 0.329411774873734],...
    'DisplayName','Mixed Silicates','EdgeColor',[0 0 0]);
set(bar1(3),'FaceColor',[1 1 0],'DisplayName','Ca Rich','EdgeColor',[0 0 0]);
set(bar1(4),'FaceColor',[0 0.498039215803146 0],...
    'DisplayName','Fresh Chlorides','EdgeColor',[0 0 0]);
set(bar1(5),'FaceColor',[0.87058824300766 0.490196079015732 0],...
    'DisplayName','Mixed Chlorides','EdgeColor',[0 0 0]);
set(bar1(6),'DisplayName','Sulphate','EdgeColor',[0 0 0]);
set(bar1(7),'FaceColor',[0.854901969432831 0.701960802078247 1],...
    'DisplayName','Gypsum','EdgeColor',[0 0 0]);
set(bar1(8),'FaceColor',[0.96078431372549 0.92156862745098 0.92156862745098],'DisplayName','Carbonaceous','EdgeColor',[0 0 0]);
set(bar1(9),'FaceColor',[1 0 0],'DisplayName','Biogenic','EdgeColor',[0 0 0]);
set(bar1(10),'FaceColor',[0 0 0],'DisplayName','Phosphate','EdgeColor',[0 0 0]);
set(bar1(11),'FaceColor',[1 0.600000023841858 0.7843137383461],...
    'DisplayName','Metallic','EdgeColor',[0 0 0]);
set(bar1(12),...
    'FaceColor',[0.494117647409439 0.494117647409439 0.494117647409439],...
    'DisplayName','Biomass Tracers','EdgeColor',[0 0 0]);
set(bar1(13),...
    'FaceColor',[0.952941179275513 0.87058824300766 0.733333349227905],...
    'DisplayName','Other','EdgeColor',[0 0 0]);
legend('Silicates','Mixed Silicates','Ca Rich','Fresh Chlorides','Mixed Chlorides',...
    'Sulphate','Gypsum','Carbonaceous','Biogenic','Phosphate',...
    'Metallic','Biomass Tracers','Other',-1) 
xlabel('Size Range (\mum)','fontsize',16)
ylabel('Fraction','fontsize',16)
set(gca,'ylim',[0 1])
titlelab=num2str(length(OUTPUT.RAW.C));
title(['n= ',titlelab]);
        
        elseif flag==2,
            
        %
%                  10 Classes in Total
%                       Other, Secondary, Sulphates, Carbonates,
%                       Phosphates, Chlorides, Oxides, Quartz, Silicates,
%                       Mixtures
%
%


       for i=1:length(lim_upper(1:end-1)),
            eval(['dust.ind',num2str(i,'%d'),'=find(OUTPUT.Silicate.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Silicate.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['quartz.ind',num2str(i,'%d'),'=find(OUTPUT.Quartz.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Quartz.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['carb.ind',num2str(i,'%d'),'=find(OUTPUT.Carbonates.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Carbonates.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['cl.ind',num2str(i,'%d'),'=find(OUTPUT.Chlorides.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Chlorides.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['oxides.ind',num2str(i,'%d'),'=find(OUTPUT.Oxides.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Oxides.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])            
            eval(['sulphate.ind',num2str(i,'%d'),'=find(OUTPUT.Sulphates.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Sulphates.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['phosphate.ind',num2str(i,'%d'),'=find(OUTPUT.Phosphates.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Phosphates.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['mix.ind',num2str(i,'%d'),'=find(OUTPUT.Mixtures.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Mixtures.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['secondary.ind',num2str(i,'%d'),'=find(OUTPUT.Secondary.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Secondary.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['other.ind',num2str(i,'%d'),'=find(OUTPUT.Other.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.Other.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])            
            eval(['totl.ind',num2str(i,'%d'),'=find(OUTPUT.All.Raw(:,3)>lim_lower(',num2str(i,'%d'),')'...
                ' & OUTPUT.All.Raw(:,3)<=lim_upper(',num2str(i,'%d'),'));'])
            eval(['ns',num2str(i,'%d'),'_bins=length(totl.ind',num2str(i,'%d'),');'])
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['inpdust_bins(',num2str(i,'%d'),')=length(dust.ind',num2str(i,'%d'),');']);
            eval(['inpquartz_bins(',num2str(i,'%d'),')=length(quartz.ind',num2str(i,'%d'),');']);
            eval(['inpcarb_bins(',num2str(i,'%d'),')=length(carb.ind',num2str(i,'%d'),');']);
            eval(['inpcl_bins(',num2str(i,'%d'),')=length(cl.ind',num2str(i,'%d'),');']);     
            eval(['inpox_bins(',num2str(i,'%d'),')=length(oxides.ind',num2str(i,'%d'),');']);         
            eval(['inpsulphate_bins(',num2str(i,'%d'),')=length(sulphate.ind',num2str(i,'%d'),');']);
            eval(['inpphos_bins(',num2str(i,'%d'),')=length(phosphate.ind',num2str(i,'%d'),');']);
            eval(['inpmix_bins(',num2str(i,'%d'),')=length(mix.ind',num2str(i,'%d'),');']);
            eval(['inpsecondary_bins(',num2str(i,'%d'),')=length(secondary.ind',num2str(i,'%d'),');']);
            eval(['inpother_bins(',num2str(i,'%d'),')=length(other.ind',num2str(i,'%d'),');']);
        end
        
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['inp',num2str(i,'%d'),'_bins=[inpdust_bins(',num2str(i,'%d'),') inpquartz_bins(',num2str(i,'%d'),') inpcarb_bins(',num2str(i,'%d'),...
                ') inpcl_bins(',num2str(i,'%d'),') inpox_bins(',num2str(i,'%d'),') inpsulphate_bins(',num2str(i,'%d'),...
                ') inpphos_bins(',num2str(i,'%d'),') inpmix_bins(',num2str(i,'%d'),') inpsecondary_bins(',num2str(i,'%d'),...
                ') inpother_bins(',num2str(i,'%d'),')];']);
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['input',num2str(i,'%d'),'_bins=inp',num2str(i,'%d'),'_bins./sum(inp',...
                num2str(i,'%d'),'_bins);']);
        end
        
        for i=1:length(lim_upper(1:end-1)),
            eval(['silicateindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(1);']);
            eval(['quartzindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(2);']);
            eval(['carbindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(3);']);
            eval(['clindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(4);']);
            eval(['oxindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(5);']);
            eval(['sulphindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(6);']);
            eval(['phosindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(7);']);
            eval(['mixindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(8);']);
            eval(['secdryindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(9);']);
            eval(['otherindex(',num2str(i,'%d'),')=input',num2str(i,'%d'),'_bins(10);']);
        end




    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  


zz=figure();
bar1=bar(1:length(lim_lower(1:end-1)),[silicateindex' quartzindex' carbindex' clindex' oxindex'...
    sulphindex' phosindex' mixindex' secdryindex' otherindex'],1.0,'stack');
if case_label==1, set(gca,'xtick',1:4:length(lim_lower)); end
if case_label==2; set(gca,'xtick',1:2:length(lim_lower)); end
set(gca,'XTickLabel',ll,'fontsize',16) 
set(bar1(1),'FaceColor',[0 0 1],'DisplayName','Silicates','EdgeColor',[0 0 0]);
set(bar1(2),...
    'FaceColor',[0.34901961684227 0.200000002980232 0.329411774873734],...
    'DisplayName','Quartz','EdgeColor',[0 0 0]);
set(bar1(3),'FaceColor',[1 1 0],'DisplayName','Carbonates','EdgeColor',[0 0 0]);
set(bar1(4),'FaceColor',[0 0.498039215803146 0],...
    'DisplayName','Chlorides','EdgeColor',[0 0 0]);
set(bar1(5),'FaceColor',[1 0.600000023841858 0.7843137383461],...
    'DisplayName','Oxides','EdgeColor',[0 0 0]);
set(bar1(6),'FaceColor',[0 1 1],'DisplayName','Sulphates','EdgeColor',[0 0 0]);
set(bar1(7),'FaceColor',[0 0 0],'DisplayName','Phosphates','EdgeColor',[0 0 0]);
set(bar1(8),'FaceColor',[1 0.694117647058824 0.392156862745098],...
    'DisplayName','Mixtures','EdgeColor',[0 0 0]);
set(bar1(9),'FaceColor',[1 0 0],'DisplayName','Secondary','EdgeColor',[0 0 0]);
set(bar1(10),...
    'FaceColor',[0.992156862745098 0.917647058823529 0.796078431372549],...
    'DisplayName','Other','EdgeColor',[0 0 0]);
legend('Silicates','Quartz','Carbonates','Chlorides','Oxides',...
    'Sulphates','Phosphates','Mixtures','Secondary','Other',-1) 
xlabel('Size Range (\mum)','fontsize',16)
ylabel('Fraction','fontsize',16)
set(gca,'ylim',[0 1])
titlelab=num2str(length(OUTPUT.RAW.C));
title(['n= ',titlelab]);
            
        end

        if flag==1,
            OUTPUT.SIZESEG.CLASSIFICATIONS=[mindustindex' mixindex' carichindex' freshssindex' agedssindex' sulphateindex' gypsumindex' ...
                carbonindex' bioindex' phosindex' oxindex' biomassindex' otherindex'];
            OUTPUT.SIZESEG.BIN_UPPER=lim_upper;
            OUTPUT.SIZESEG.BIN_LOWER=lim_lower;
            OUTPUT.SIZESEG.BIN_CENT=ll;            
            OUTPUT.SIZESEG.CATEGORIES={'Silicate','MixedSil','CaRich','FreshCl','AgedCl',...
                'Sulph','Gypsum','Carbon','Bio','Phos',...
                'Metal','Biomass','Other'};
        elseif flag==2,
            OUTPUT.SIZESEG.CLASSIFICATIONS=[silicateindex' quartzindex' carbindex' clindex' oxindex'...
                sulphindex' phosindex' mixindex' secdryindex' otherindex'];
            OUTPUT.SIZESEG.BIN_UPPER=lim_upper;
            OUTPUT.SIZESEG.BIN_LOWER=lim_lower;   
            OUTPUT.SIZESEG.BIN_CENT=ll;
            OUTPUT.SIZESEG.CATEGORIES={'Silicate','Quartz','Carbonates','Chlorides','Oxides',...
                'Sulphates','Phosphates','Mixtures','Secondary','Other'};
        end
        
end
