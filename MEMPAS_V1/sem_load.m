% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % %           
% % % % % % %               LOAD IN SEM DATA
% % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

dir=pwd;
filename = uigetfile({'*.*','All Files (*.*)'},'Select SEM Output File (.csv for old ESEM, .xlsx for new ESEM)');

% % % % 
% % % %     CHOOSE INSTRUMENT
% % % % 


instrument_switch = input('Old SEM [XL30] (1), or new SEM (2):  ');
switch instrument_switch
    case 1
        disp('OLD SEM [XL30] SELECTED')
        new=false;
        instrument_string={'OLD SEM [XL30] SELECTED'};
    case 2
        disp('NEW SEM SELECTED')
        new=true;
        instrument_string={'NEW SEM SELECTED'};
end


% % % % 
% % % %     LOAD IN DATA
% % % % 

    
if new==0,
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % %
    % % % % % % %              CHOOSE ELEMENT LIST
    % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    disp(['For Old SEM [XL30], element list must be specified'])
    disp(['Choose 1 of the following:'])
    disp('Young et al., 2016: [C O Na Mg Al Si P S Cl K Ca Ti Cr Fe Ni Cu Zn]')
    disp('Kandler et al., 2011: [C O Na Mg Al Si P S Cl K Ca Ti Cr Mn Fe]')
    
    [num]=csvread(filename,13,0);
    element_list=input('Choose element list: (1) Young et al., 2016, or (2) Kandler et al., 2011: ');
    switch element_list
        case 1
            disp('[C O Na Mg Al Si P S Cl K Ca Ti Cr Fe Ni Cu Zn]')
            txt={'Part#','Field#','Phase#','X_stage','Y_stage','X_cent','Y_cent','X_left','Y_low','X_width','Y_height',...
                'Xferet','Yferet','AvgDiam','LProj','Area','Perim','Shape','Aspect Ratio','Orientation','C','O','Na',...
                'Mg','Al','Si','P','S','Cl','K','Ca','Ti','Cr','Fe','Ni','Cu','Zn','CPS','AvgVideo'};
            element_string1={'Elements chosen:'};
            element_string2={'[C O Na Mg Al Si P S Cl K Ca Ti Cr Fe Ni Cu Zn]'};
        case 2
            disp('[C O Na Mg Al Si P S Cl K Ca Ti Cr Mn Fe]')
            txt={'Part#','Field#','Phase#','X_stage','Y_stage','X_cent','Y_cent','X_left','Y_low','X_width','Y_height',...
                'Xferet','Yferet','AvgDiam','LProj','Area','Perim','Shape','Aspect Ratio','Orientation','C','O','Na',...
                'Mg','Al','Si','P','S','Cl','K','Ca','Ti','Cr','Mn','Fe','CPS','AvgVideo'};
            element_string1={'Elements chosen: '};
            element_string2={'[C O Na Mg Al Si P S Cl K Ca Ti Cr Mn Fe]'};
        otherwise
            disp('OPTION NOT AVAILABLE')
    end
    x=size(num);col=x(2);leng=x(1);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % %
    % % % % % % %              ORGANISE SEM DATA
    % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    last_col=7+((col-2)-21);
    data=zeros(length(num(:,1)),last_col);
    headers=[txt(2) txt(16) txt(14) 'Volume' txt(18) txt(19) txt(21:col-2)];
    data(:,1)=num(:,2);                      % FIELD NUMBER
    data(:,2)=num(:,16);                     % PARTICLE AREA (UM^2)
    data(:,3)=num(:,14);                     % PARTICLE AVERAGE DIAMETER (UM)
    data(:,4)=pi./6.*(data(:,2).^3);         % VOLUME OF EQUIVALENT-DIAMETER SPHERE
    data(:,5)=num(:,18);                     % SHAPE PARAMETER
    data(:,6)=num(:,19);                     % ASPECT RATIO
    data(:,7:last_col)=num(:,21:col-2);      % ELEMENT LIST
    
    %     clear col leng
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % %
    % % % % % % %              APPLY CLASSIFICATION SCHEME
    % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    classifications = input('Choose classification scheme: (1) Young et al., 2016, (2) Kandler et al., 2011:  ');
    switch classifications
        case 1
            artefacts_switch1 = input('Use artefact removal?: Yes (1), or No (2): ');
            switch artefacts_switch1
                
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % %
                % % % % % % %              ARTEFACT REMOVAL
                % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                
                case 1
                    
                    disp('YOUNG ET AL., (2016) SELECTED WITH ARTEFACT REMOVAL')
                    [dat,ARTEFACTS,total_num]=sem_artefacts(data,headers);
                    [DATA]=sem_classifications(dat,headers,last_col);
                    
                    disp(['Total Number of Particles Detected: ',num2str(length(num(:,1)))])
                    disp(['Number of Viable Particles: ',num2str(length(dat(:,1)))])
                    
                    flag=1;
                    artefact_string={'YOUNG ET AL., (2016) SELECTED WITH ARTEFACT REMOVAL'};
                    
                case 2
                    
                    disp('YOUNG ET AL., (2016) SELECTED, NO ARTEFACT REMOVAL')
                    [DATA]=sem_classifications(data,headers,last_col);
                    
                    disp(['Total Number of Particles Detected: ',num2str(length(num(:,1)))])
                    disp(['Number of Viable Particles: ',num2str(length(data(:,1)))])
                    
                    flag=1;
                    artefact_string={'YOUNG ET AL., (2016) SELECTED, NO ARTEFACT REMOVAL'};
                    
                otherwise
                    disp('OPTION NOT AVAILABLE')
                    
                    
            end
        case 2
            artefacts_switch2 = input('Use artefact removal?: Yes (1), or No (2): ');
            switch artefacts_switch2
                
                case 1
                    
                    disp('KANDLER ET AL., (2011) SELECTED WITH ARTEFACT REMOVAL')
                    [dat,ARTEFACTS,total_num]=sem_artefacts(data,headers);
                    [DATA]=sem_kandler2011(dat,headers,last_col);
                    
                    disp(['Total Number of Particles Detected: ',num2str(length(num(:,1)))])
                    disp(['Number of Viable Particles: ',num2str(length(dat(:,1)))])
                    
                    flag=2;
                    artefact_string={'KANDLER ET AL., (2011) SELECTED WITH ARTEFACT REMOVAL'};
                    
                case 2
                    
                    disp('KANDLER ET AL., (2011) SELECTED, NO ARTEFACT REMOVAL')
                    [DATA]=sem_kandler2011(data,headers,last_col);
                    
                    disp(['Total Number of Particles Detected: ',num2str(length(num(:,1)))])
                    disp(['Number of Viable Particles: ',num2str(length(data(:,1)))])
                    
                    flag=2;
                    artefact_string={'KANDLER ET AL., (2011) SELECTED, NO ARTEFACT REMOVAL'};
                    
                otherwise
                    disp('OPTION NOT AVAILABLE')
            end
            
    end
    
    
    
elseif new==1
    
% % %     [~,leng] = xlsread(filename,'Mass % (norm.)', 'A:A'); %number of rows in array
% % %     [~,col] = xlsread(filename,'Mass % (norm.)', '7:7'); %number of columns in array
    [~,leng] = xlsread(filename,5, 'A:A'); %number of rows in array
    [~,col] = xlsread(filename,5, '7:7'); %number of columns in array    
%     Alphabet1=('A':'Z').';
% % %     for i=1:26, Alphabet2{i}=strcat('A',Alphabet1(i)); end;
%     a=repmat(char('A'),length(Alphabet1),1);
%     Alphabet2=strcat(a,Alphabet1);    
%     if length(col)>26 && length(col)<52,
%         col_num=Alphabet2(length(col)-26,:);
%     else col_num=Alphabet1(length(col));
%     end
%     end_col=[col_num,num2str(length(leng))];
    
% %     [num,txt,raw]=xlsread(filename,'Mass % (norm.)',['B7:',end_col]);
% %     [num,txt,raw]=xlsread(filename,5,['B7:',end_col]);

%     lastrow=num2str(length(leng));
%     [field_no]=xlsread(filename,5,['B7:B',lastrow]);
%     [area]=xlsread(filename,5,['F7:F',lastrow]);
%     [avg_diam]=xlsread(filename,5,['G7:G',lastrow]);
%     [volume]=xlsread(filename,5,['L7:L',lastrow]);
    

    filename_csv=[filename(1:end-5),'.csv'];
    element_list=col(16:end);
    [num]=csvread(filename_csv,7,15);       % read in elements
    txt=[element_list 'Field No' 'Area' 'Average Diameter' 'Volume'];
%     end_row=perl('countlines.pl',filename);
%     [num]=csvread(filename,8,15);

   
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % %
    % % % % % % %              ORGANISE SEM DATA
    % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
    
% %     last_col=5+(length(col)-16);
% %     data=zeros(length(num(:,1)),last_col);
% %     headers=[txt(1,1) txt(1,5) txt(1,6) txt(1,11) txt(1,15:length(col)-1)];
% %     data(:,1)=num(:,1);
% %     data(:,2:3)=num(:,5:6);
% %     data(:,4)=num(:,11);
% %     data(:,5:last_col)=num(:,15:length(col)-1);
% %     data(isnan(data))=0;
    

    last_col=5+(length(col)-16);
    data=zeros(length(num(:,1)),last_col);
    headers=[txt(end-3) txt(end-2) txt(end-1) txt(end) txt(1:end-4)];
    data(:,1)=num(:,end-3);
    data(:,2)=num(:,end-2);
    data(:,3)=num(:,end-1);
    data(:,4)=num(:,end);
    data(:,5:last_col)=num(:,1:last_col-4);
    data(isnan(data))=0;


    element_string1={'Elements identified: '};
% %     element_string2=headers(5:last_col);
    element_string2=element_list;
    
    disp(['Total Number of Particles Detected: ',num2str(length(num(:,1)))])
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % %
    % % % % % % %              APPLY CLASSIFICATION SCHEME
    % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    classifications = input('Choose classification scheme: (1) Young et al., 2016, (2) Kandler et al., 2011:  ');
    switch classifications
        case 1
            disp('YOUNG ET AL., (2016) SELECTED')
            flag=1;
            [DATA]=sem_classifications(data,headers,last_col);
            artefact_string={'YOUNG ET AL., (2016) SELECTED'};
            
        case 2
            disp('KANDLER ET AL., (2011) SELECTED')
            [DATA]=sem_kandler2011(data,headers,last_col);
            flag=2;
            artefact_string={'KANDLER ET AL., (2011) SELECTED'};
            
        otherwise
            disp('OPTION NOT AVAILABLE')
    end
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % %
% % % % % % % % % % %     SIZE SEGREGATED PARTICLE COMPOSITIONS
% % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% fig_switch = input('Output size-segregated particle composition figure?: Yes (1), No (2): ');
% switch fig_switch
%     case 1
disp('SIZE SEGREGATION:')
[DATA,fig_handle1]=sem_sizeseg(DATA,flag);
fig_string={'SIZE SEGREGATION COMPLETED'};
%     case 2
%         disp('SIZE SEGREGATION NOT SELECTED')
%         fig_string={'SIZE SEGREGATION NOT SELECTED'};
%     otherwise
%         disp('OPTION NOT AVAILABLE')
% end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     SIZE DISTRIBUTIONS
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

sizedist_switch = input('Output size distribution?: Yes (1), No (2): ');
switch sizedist_switch
    case 1
        disp('SIZE DISTRIBUTION SELECTED')
        DATA.INITIAL.VOLUME_L=input('Enter volume of air sampled by filter (/L): ');
        DATA.INITIAL.AREA_mm=input('Enter area covered by scan (mm): ');
        DATA.INITIAL.FILTERDIAM_mm=input('Enter filter size (mm): ');
        DATA.INITIAL.EXPOSURE_s=input('Enter length of exposure (s): ');
        [DATA,fig_handle2]=sem_sizedist(DATA,flag);
        sizedist_string1={'SIZE DISTRIBUTION SELECTED'};
        sizedist_flag=1;
        disp('PLOT= plot(DATA.SIZEDIST.SIZEBINS_AV,DATA.SIZEDIST.DNDLOGD_AV_cm3):')
    case 2
        disp('SIZE DISTRIBUTION NOT SELECTED')
        sizedist_string1={'SIZE DISTRIBUTION NOT SELECTED'};
        sizedist_string2={' '};
        sizedist_flag=0;
    otherwise
        disp('OPTION NOT AVAILABLE')
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     OUTPUT DATA
% % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

dirname = input('A new directory will be made for output data. Enter directory name: ','s');
mkdir(dirname)
cd(dirname)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     ORGANISE DATA FOR .CSV OUTPUT
% % % % % % % % % % % %         ELEMENTAL DATA
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if new==0
    output_data=data;
    output_headers=[headers 'Classification Flag'];
    for i=1:length(DATA.SIZESEG.CATEGORIES),
        eval(['output_data(DATA.',char(DATA.SIZESEG.CATEGORIES(i)),'.All_Indices,24)=',num2str(i),';'])
    end
    
    disp(['Enter Classification Flag (col 24) into DATA.SIZESEG.CATEGORIES for full description, or see ***_info.txt file for details'])
    
    
    % %// Write fields to CSV file
    name=input('Enter output filename for elemental data: ','s');
    outputfilename=strcat(char(name),'.csv');
    fid = fopen(outputfilename, 'w');
    fmt_str = repmat('%s,', 1, size(output_headers));
    fprintf(fid, [fmt_str(1:end-1),'\n'],output_headers{:});
    fclose(fid);
    dlmwrite(outputfilename, output_data, '-append', 'precision', '%.6f', 'delimiter', '\t');
    
elseif new==1
    output_data=data;
    output_headers=[headers 'Classification Flag'];
    for i=1:length(DATA.SIZESEG.CATEGORIES),
        eval(['output_data(DATA.',char(DATA.SIZESEG.CATEGORIES(i)),'.All_Indices,last_col+1)=',num2str(i),';'])
    end
    
    disp(['Enter Classification Flag (col 24) into DATA.SIZESEG.CATEGORIES for full description, or see ***_info.txt file for details'])
    
    
    % %// Write fields to CSV file
    name=input('Enter output filename for elemental data: ','s');
    outputfilename=strcat(char(name),'.csv');
    fid = fopen(outputfilename, 'w');
    fmt_str = repmat('%s,', 1, size(output_headers));
    fprintf(fid, [fmt_str(1:end-1),'\n'],output_headers{:});
    fclose(fid);
    dlmwrite(outputfilename, output_data, '-append', 'precision', '%.6f', 'delimiter', '\t');
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     ORGANISE DATA FOR .CSV OUTPUT
% % % % % % % % % % % %             SIZE DISTRIBUTIONS
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if sizedist_flag==1,
    for i=1:length(DATA.SIZESEG.CATEGORIES)
        eval(['zeroflag(i)=length(DATA.',char(DATA.SIZESEG.CATEGORIES(i)),'.All_Indices);'])
    end
    index1=find(zeroflag==0);
    for i=1:length(index1),
        eval(['DATA.',char(DATA.SIZESEG.CATEGORIES(index1(i))),'.NumberConcentration_cm3=zeros(length(DATA.SIZEDIST.SIZEBINS),1)';])
    end
    index2=find(zeroflag==1);
    for i=1:length(index2),
        eval(['DATA.',char(DATA.SIZESEG.CATEGORIES(index2(i))),'.NumberConcentration_cm3=zeros(length(DATA.SIZEDIST.SIZEBINS),1)';])
    end
    
%     sizedist_switch2 = input('Choose data to output to .csv file: (1) Number concentration/size bin, (2) dN/dlogDp: ');
%     switch sizedist_switch2
%         case 1
            sizedist_headers=['Bin centres (um)','Total_cm3',DATA.SIZESEG.CATEGORIES];
            sizedist_data(:,1)=DATA.SIZEDIST.SIZEBINS;
            sizedist_data(:,2)=DATA.SIZEDIST.NUMCONC_cm3;
            for i=1:length(DATA.SIZESEG.CATEGORIES)
                eval(['sizedist_data(:,2+i)=DATA.',char(DATA.SIZESEG.CATEGORIES(i)),'.NumberConcentration_cm3;'])
            end
            sizedist_string2={'OUTPUT: NUMBER CONCENTRATION (/CM3) IN EACH SIZE BIN'};
%         case 2
%             sizedist_headers=['Bin centres (um)','Total_cm3',DATA.SIZESEG.CATEGORIES];
%             sizedist_data(:,1)=DATA.SIZEDIST.SIZEBINS(1:end-1);
%             sizedist_data(:,2)=DATA.SIZEDIST.DNDLOGD_cm3;
%             for i=1:length(DATA.SIZESEG.CATEGORIES)
%                 eval(['sizedist_data(:,2+i)=DATA.',char(DATA.SIZESEG.CATEGORIES(i)),'.DNDLOGD_cm3;'])
%             end 
%             sizedist_string2={'OUTPUT: DN/DLOGDP'};
%     end

    
    
    % %// Write fields to CSV file
    sizedistname=input('Enter output filename for size distribution data: ','s');
    sizedistfilename=strcat(sizedistname,'.csv');
    fid = fopen(sizedistfilename, 'w');
    fmt_str = repmat('%s,', 1, size(sizedist_headers));
    fprintf(fid, [fmt_str(1:end-1),'\n'],sizedist_headers{:});
    fclose(fid);
    dlmwrite(sizedistfilename, sizedist_data, '-append', 'precision', '%.6f', 'delimiter', '\t');
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     ORGANISE INFO FOR .TXT FILE
% % % % % % % % % % % %             METADATA
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% %// Write metadata to text file
string={'Classification flag corresponds to the following:'};
string_blank={' '};
for i=1:length(DATA.SIZESEG.CATEGORIES),
   string1(i)=strcat(num2str(i),': ',DATA.SIZESEG.CATEGORIES(i));
end

metafilename=strcat(name,'_info.txt');
fid = fopen(metafilename, 'w');
fmt_str = [repmat('%s\n', 1, 50),'%s\n'];
fprintf(fid, fmt_str,instrument_string{:});
fprintf(fid, fmt_str,element_string1{:});
fprintf(fid, fmt_str,element_string2{:});
fprintf(fid, fmt_str,artefact_string{:});
fprintf(fid, fmt_str,fig_string{:});
fprintf(fid, fmt_str,sizedist_string1{:});
fprintf(fid, fmt_str,sizedist_string2{:});
fprintf(fid, fmt_str,string_blank{:});
fprintf(fid, fmt_str,string{:});
fprintf(fid, fmt_str,string_blank{:});
% fprintf(fid, [fmt_str(1:end-1),'\n'],string{:});
for i=1:length(DATA.SIZESEG.CATEGORIES),
    fprintf(fid,[fmt_str(1:end-1),'\n'],string1{i});
end
fclose(fid);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % %
% % % % % % % % % % % %     CLEAR WORKSPACE
% % % % % % % % % % % %     SAVE MATLAB STRUCTURE
% % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clearvars -except DATA
uisave

cd ..
