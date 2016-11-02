%%
%
%       ARTEFACT REMOVAL FUNCTION
%           USED WITH OLD SEM DATA ONLY
%

function [good_num,ART,total_num]=sem_artefacts(num,txt)

particle_size=num(:,3);          % sets size measurements to a variable
total_num=num;
unclassed=num;

m=length(num);

%
%       Ensure all real particles
%

fix=find(particle_size >= 0.12);
good_num=unclassed(fix,:);
ind0=find(nansum(num(:,7:end)) == 0);       % MAKE SURE COMPOSITIONAL DATA IS NON-ZERO
indsh=find(num(:,5) >=30);                  % EXCLUDE ABNORMAL SHAPES
indasp=find(num(:,6) >=30);                 % EXCLUDE ABNORMAL ASPECT RATIOS


good_num=unclassed;
n=length(good_num);

total=zeros(n,1);
for i=1:n,
    comp=good_num(i,:);
    total(i)=sum(good_num(i,9:end));    % SUM OF ELEMENTS >= NA
end

% % % % % % % % % % % % %
% % % % % % % % % % % % %
% % % % % % % % % % % % %               Remove Filter Artefacts
% % % % % % % % % % % % %
% % % % % % % % % % % % %
artefact1=zeros(n,length(txt));
artefact2=zeros(n,length(txt));
artefact3=zeros(n,length(txt));
artefact4=zeros(n,length(txt));
artefact5=zeros(n,length(txt));
artefact6=zeros(n,length(txt));
flag_zero=zeros(n,length(good_num(1,9:end)));
flag_lt05=zeros(n,length(good_num(1,9:end)));
flag_lt15=zeros(n,length(good_num(1,9:end)));

       
        
for i=1:n,
    for j=1:length(good_num(1,9:end)),
        if good_num(i,8+j) == 0;    % equal to 0 wt%
            flag_zero(i,j)=1;
        else flag_zero(i,j)=0;
        end
        if good_num(i,8+j)<=0.5;    % less than 0.5 wt%
            flag_lt05(i,j)=1;
        else flag_lt05(i,j)=0;
        end
        if good_num(i,8+j)<=1.2;    % less than 1.5 wt%
            flag_lt15(i,j)=1;
        else flag_lt15(i,j)=0;
        end
    end
end

if length(find(strcmp('C', txt)==1))==1, C=num(:,find(strcmp('C', txt)==1)); C(isnan(C))=0; end
if length(find(strcmp('O', txt)==1))==1, O=num(:,find(strcmp('O', txt)==1)); O(isnan(O))=0; end
if length(find(strcmp('Na', txt)==1))==1, Na=num(:,find(strcmp('Na', txt)==1)); Na(isnan(Na))=0; Na_index=find(strcmp('Na', txt)==1); else Na=zeros(n,1); Na_index=0; end
if length(find(strcmp('Fe', txt)==1))==1, Fe=num(:,find(strcmp('Fe', txt)==1)); Fe(isnan(Fe))=0; Fe_index=find(strcmp('Fe', txt)==1); else Fe=zeros(n,1); Fe_index=0; end
if length(find(strcmp('Ni', txt)==1))==1, Ni=num(:,find(strcmp('Ni', txt)==1)); Ni(isnan(Ni))=0; Ni_index=find(strcmp('Ni', txt)==1); else Ni=zeros(n,1); Ni_index=0; end
if length(find(strcmp('Cu', txt)==1))==1, Cu=num(:,find(strcmp('Cu', txt)==1)); Cu(isnan(Cu))=0; Cu_index=find(strcmp('Cu', txt)==1); else Cu=zeros(n,1); Cu_index=0; end
if length(find(strcmp('Zn', txt)==1))==1, Zn=num(:,find(strcmp('Zn', txt)==1)); Zn(isnan(Zn))=0; Zn_index=find(strcmp('Zn', txt)==1); else Zn=zeros(n,1); Zn_index=0; end

for i=1:n,
    %   Zn>Cu>Ni>Fe and C>=70 and Size<0.3
    if Zn(i) > Cu(i)
        if Cu(i) > Ni(i)
            if Ni(i) > Fe(i)
                if good_num(i,3) <= 0.3
                    if C(i) >= 70
                        if sum(flag_zero(i,:),2)<=4
                            artefact1(i,:)=good_num(i,:);
                        end
                    end
                end
            end
        end
    end
    
    thresh=0.97*nansum(num(i,7:end));
    %   C+O>97%
    if C(i) + O(i) > thresh
        if sum(flag_zero(i,:),2)<=6         % less than 6 are equal to 0
            artefact2(i,:)=good_num(i,:);
        end
    end
    
    
    %
    %       If none are equal to zero:
    %
    flag0=sum(flag_zero(i,:),2);    %       If none are equal to zero
    if flag0 == 0;
        if C(i) >= 70
            if good_num(i,3) <= 0.3
                artefact3(i,:)=good_num(i,:);
            end
        end
    end
    clear flag0
end


%
%       If all are <2 wt% except Ni, Cu and Zn which may be zero
%

for i=1:n,
    flag0=sum(flag_zero(i,:),2);                    % sum of number of zeros in row
    
    if Fe_index>0,                                  % if Fe is present
        flag1=sum(flag_zero(i,1:Fe_index-8),2);         % flag1 = sum number of zeros in Na:Fe 
    elseif size(flag_zero,2)<Fe_index-8             % if Fe is not present,      
        flag1=sum(flag_zero(i,1:end),2);                % flag1 = sum number of zeros to heaviest element
    end 
    
        if Fe_index<size(good_num,2)
            if good_num(i,9:end)<=2,
                if good_num(i,9:end)>=0 
                    artefact4(i,:)=good_num(i,:);
                elseif flag1<=3  
                artefact4(i,:)=good_num(i,:);
                end
            end
        elseif Fe==size(good_num,2)
            artefact4(i,:)=0;
        end        
        
    %
    %           If all are <0.5 %wt or <1.0 wt%
    %
    flag2=sum(flag_lt05(i,:),2);     % less than 0.5 wt%
    flag3=sum(flag_lt15(i,:),2);       % less than 1.5 wt%
    if flag2 == length(good_num(1,9:end));
        if flag0<=5
            if C(i) >= 60
                artefact5(i,:)=good_num(i,:);
            end
        end
    end
    if flag3 == length(good_num(1,9:end));
        if C(i) >= 60
            if flag0<=6
                artefact6(i,:)=good_num(i,:);
            end
        end
    end
    clear flag0 flag1 flag2 flag3
end

artef1_index=find(artefact1(:,3)>0);
artef2_index=find(artefact2(:,3)>0);
artef3_index=find(artefact3(:,3)>0);
artef4_index=find(artefact4(:,3)>0);
artef5_index=find(artefact5(:,3)>0);
artef6_index=find(artefact6(:,3)>0);

hh=[artef1_index;artef2_index;artef3_index;artef4_index;artef5_index;artef6_index;ind0;indsh;indasp];
artef_index=unique(hh);
ART.TOTAL=good_num(artef_index,:);

if isempty(artef_index)
    disp('No Artefacts Detected')
else
    good_num(artef_index,:)=0;
end

good_num(artef_index,:)=[];

end

