%% Indexing Module (DRM) Demo
%    Copyright (C) 2012  Juan M. Banda, Rafal A. Angryk from Montana State University
%    Contact: juan@jmbanda.com
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% MORE INFO:
%  Stand-alone script to ilustrate the usage of the Feature Extraction
%  Module. Remember to set your path to the correct place where the demo
%  subset of the dataset has been extracted
%  For more details
%  Juan M. Banda's dissertation:
%  "FRAMEWORK FOR CREATING LARGE-SCALE CONTENT-BASED IMAGE RETRIEVAL SYSTEM
%  (CBIR) FOR SOLAR DATA ANALYSIS"
%  http://www.jmbanda.com/Dissertation/
%
%  Notes on this DEMO:
%  http://www.jmbanda.com/Framework/Demo/
%
%  Note: This module will generate the necesary files to be used with the
%  modified implementation of GIMP in C (to be used on a Linux evironment

clearvars;
clases=8;                           %Number of classes
Csize=20;                            %Number of images in each folder
total_im=clases*Csize;
segments=8;                         %Number of Columns/rows  NbyN
class_n=20;                         %class size
sample_size=clases*class_n;         %TOTAL number of images
cells=segments*segments;            %Grid Cell Size
parameters=10;                      %Number of Image Parameters
distances=12;
norm=1;
index=segments*segments*parameters; %Number of elements in the Feature Vector
dataSet='Solar';                    %Dataset Being Manipulated
pathh='C:\Demo\';                   %Root folder of datasets
%Save too path
pathhTOS='C:\Demo\';                   %Root folder of datasets
exten='tif';
clasesNames={'Active Region' 'Coronal Jet' 'Emerging Flux' 'Filament' 'Filament Activation' 'Filament Eruption' 'Flare' 'Oscillation'};
exlen=length(char(exten))+1;
plt=0;       %0 for plots , 1 for no plots
% IF WE WANT ALL dimensions
ALL_DIM=1;   %1 yes, 0 no
% IF WE WANT DIMENSIONALLY REDUCED
DIM_RED=0;   %1 yes, 0 no
% IF WE WANT TO RUN Dimensionality reduction in thesis context
DIM_RED_T=1;
%Image Parameter List 
imgParam={'Entropy' 'Mean' 'StandardDev' 'FracDim' 'Moment3' 'Moment4' 'Uniformity' 'RelSm' 'TamDir' 'TamCont'};
%% This is form the DRM Module, in the demo we hardcoded the results
dim_to=[15,18,22,28,31,39,51,72];
n_DM=size(dim_to,2); %Number of dimensions
%%%%%%%%%%%%%%%%%%%%%% END OF GLOBAL VARIABLES
%% Get Data from the FE files to manipulate
Histogram_Array=zeros(parameters,Csize*clases,segments*segments);  %Pre-allocate
imageCounter=1;
clear labelsData;
for fv=1:clases
folder=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\*',exten)]);
folder2=char([strcat(pathh,dataSet,'\',clasesNames(fv),'\')]);
list=dir(folder);
    for ii=1:class_n
        filenam=strcat(folder2,list(ii,1).name);
        s2=length(char(filenam));
        filenam=strcat(filenam(1:s2-exlen),'-',int2str(segments),'x',int2str(segments),'.txt');
        [sample2(:,1),sample2(:,2),sample2(:,3),sample2(:,4),sample2(:,5),sample2(:,6),sample2(:,7),sample2(:,8),sample2(:,9),sample2(:,10),sample2(:,11)] = textread(filenam,'%s %s %s %s %s %s %s %s %s %s %s','delimiter', '\t','bufsize', 8191);
        for sgems=1:cells
            for parmsts=1:10
                Histogram_Array(parmsts,imageCounter,sgems)=str2double(sample2(sgems,parmsts));
            end
            labelsData(imageCounter)=sample2(sgems,11);
        end
        imageCounter=imageCounter+1;
    end %Looop files in class
end %loop classes
%Switch to images , segments*segments*parameters representation
imsr=1;
imag=1;
tmpH=zeros(total_im,(cells*parameters));
for imgk=1:total_im
    for tmy=1:cells
        for cont=1:parameters
            tmpH(imsr,(((tmy-1)*10)+cont))=Histogram_Array(cont,imag,tmy);
        end
    end    
    imag=imag+1;
    imsr=imsr+1;
end

%NORMALIZE STUFF
if norm==1
    clear maxM;
    clear minM;
    maxM =zeros(parameters,1);
    minM =zeros(parameters,1);
    %Normalize  between 0 and 1
    for parm=1:parameters
        maxM(parm)=max(max(Histogram_Array(parm,:,:)));  % Max of Image Param 1
        minM(parm)=min(min(Histogram_Array(parm,:,:)));  % Max of Image Param 1
    end
    for nmH=1:total_im
        for tmy=1:cells
            for clse=1:parameters
                Histogram_Array(clse,nmH,tmy)=(Histogram_Array(clse,nmH,tmy) + abs(minM(clse))) / (maxM(clse)+abs(minM(clse)));
            end
        end  
    end
end %Normalization

%% Write ALL Dimensions (original Dataset)
if ALL_DIM==1
    %SAVE THE NORMALIZED MAPPED SPACE
    %Binary mode of Data points
    clusters=round(index/2);  %Twice the dimensionality
    %W
    fill1=strcat(dataSet,'-D-',num2str(index));
    fil1=strcat(pathhTOS,fill1,'-CLU-',num2str(clusters),'.bin');
    fid = fopen(fil1, 'w');
    for tmpy=1:total_im
        fwrite(fid, tmpH(tmpy,:), 'float');
    end
    fclose(fid);
    %FIND REFERENCE POINTS
    %Clustering Section
    try
        [IDX,C]=kmeans(tmpH,clusters);
    catch ER
        clusters=round(clusters/2); 
        try 
            [IDX,C]=kmeans(tmpH,clusters);
        catch ER2
            clusters=round(clusters/2); 
            try 
                [IDX,C]=kmeans(tmpH,clusters);
            catch ER3
                clusters=round(clusters/2); 
                [IDX,C]=kmeans(tmpH,clusters);
            end
        end
    end
    %IDX contains the cluster the points belong
    %C contains the cluster locations
    %Sequential mode of reference points
    fil2=strcat(pathhTOS,dataSet,'-D-',num2str(index),'-CLU-',num2str(clusters),'.ref');
    fid = fopen(fil2, 'w');
    for tmpy=1:clusters
        fprintf(fid,'%f ', C(tmpy,:));
        fprintf(fid,'\n');
    end
    fclose(fid);        
end
%% Splitting the dataset in 67 - 33% samples (equally balanced....) this
% might need work for different datasets
% basically here the code selected 2 instances for training and the
% third one for testing (and repeats)
% Normalize - before separating
for nmH=1:total_im
    intBy=sum(tmpH(nmH,:));
    for tmy=1:(cells*parameters)
        tmpH(nmH,tmy)=tmpH(nmH,tmy)/intBy;
        if tmpH(nmH,tmy)==0
            tmpH(nmH,tmy)=0.000000000000000000000001;
        end
    end
end   
clear tmpH1;
clear tmpH2;
imag=1;
imsr=1;
cont=1;
imSRR=1;
for imgk=1:total_im
    %Test
    tmpH1(imsr,:)=tmpH(imag,:);
    labelsData1(imsr)=labelsData(imag);
    imsr=imsr+1;
    imag=imag+1;
    imgk=imgk+1;
    if imag==total_im+1
        break
    end
    tmpH1(imsr,:)=tmpH(imag,:);
    labelsData1(imsr)=labelsData(imag);
    imsr=imsr+1;    
    imgk=imgk+1;
    imag=imag+1;
    if imag==total_im+1
        break
    end
    tmpH2(imSRR,:)=tmpH(imag,:);
    labelsData2(imSRR)=labelsData(imag);
    imSRR=imSRR+1;
    imag=imag+1;
    if imag==total_im+1
        break
    end
end
[sH1_1 sH1_2]=size(tmpH1);
[sH2_1 sH2_2]=size(tmpH2);
%% Actual section of dimensionality reduction
if DIM_RED==1
 for lp_n_d=1:n_DM  %Loop through the targeted dimensions
 no_dimsEV=dim_to(lp_n_d);
 % Loop for 6 dim red techniques
 for dmr_t=1:8
    skip=0;
    if dmr_t==1
        %PCA STUFF
        methd='PCA';
        clear mappedX;
        clear t_points;
        [mappedX, mappingPCA] = compute_mapping(tmpH1, 'PCA', no_dimsEV);
        %Generate New test set mapping
        t_points = out_of_sample(tmpH2, mappingPCA);
    end
    if dmr_t==2
        methd='SVD';
        clear U;
        clear S;
        clear V;
        [U,S,V] = svd(zscore(tmpH1),'econ');
        variances = diag(S).^2 / (size(tmpH1,1)-1);
        varExplained = 100 * variances./sum(variances);
        Variance=cumsum(variances)./sum(variances);    %base 100          
        clear mappedX;
        clear t_points;
        %Training        
        mappedX = tmpH1 * V(:,1:no_dimsEV);    
        %Test Data
        t_points= tmpH2 * V(:,1:no_dimsEV);    
    end    
    if dmr_t==3
        methd='KernelPCA';
        clear mappedX;
        clear t_points;
        [mappedX, mappingKPCA] = compute_mapping(tmpH1, methd, no_dimsEV);
        t_points = out_of_sample(tmpH2, mappingKPCA);
    end
    if dmr_t==4
        %FactorAnalysis Stuff
        methd='FactorAnalysis';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingFA] = compute_mapping(tmpH1, methd, no_dimsEV);
        %Generate New test set mapping
        t_points = out_of_sample(tmpH2, mappingFA);
    end
    if dmr_t==5
        methd='LLE';
        clear mappedX;
        clear t_points;    
        try
            [mappedX, mappingLLE] = compute_mapping(tmpH1, methd, no_dimsEV);
        catch
            msg=char(strcat('LEE-',num2str(no_dimsEV),' failed!'))
            skip=1;
        end          
        if skip ~=1
            %Generate New test set mapping
            mappedX = out_of_sample(tmpH1, mappingLLE);
            t_points = out_of_sample(tmpH2, mappingLLE);
        end
    end
    if dmr_t==6
        methd='Laplacian';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingLAP] = compute_mapping(tmpH1, methd, no_dimsEV);
        mappedX = out_of_sample(tmpH1, mappingLAP);
        t_points = out_of_sample(tmpH2, mappingLAP);
    end
    if dmr_t==7
        methd='Isomap';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingISO] = compute_mapping(tmpH1, methd, no_dimsEV);
        %Generate New test set mapping
        mappedX = out_of_sample(tmpH1, mappingISO);
        t_points = out_of_sample(tmpH2, mappingISO);
    end
    if dmr_t==8
        methd='LPP';
        clear mappedX;
        clear t_points;        
        try
            [mappedX, mappingLPP] = compute_mapping(tmpH1, methd, no_dimsEV);
        catch
            msg=char(strcat('LPP-',num2str(no_dimsEV),' failed because number of dimensions exceeds number of samples'))
            skip=1;
        end
        %Generate New test set mapping
        if skip ~=1
            t_points = out_of_sample(tmpH2, mappingLPP);
        end
    end
    if skip ~= 1
    [n1,n2]= size(mappedX);
    tst=no_dimsEV;
    if no_dimsEV ~= n2
        tst=no_dimsEV;    %Random times when something goes wrong
        no_dimsEV = n2;
    end
    NumberComp=no_dimsEV;

    %WRITE TO INDEX FILES
    %Take out the bad apples  (LLE has many!)
    smx=size(mappedX,1);
    miNew=abs(min(min(mappedX)));
    maNew=max(max(mappedX));
    normby=maNew+miNew;
    for nj=1:smx
        for nm=1:NumberComp
            mappedX(nj,nm)=(mappedX(nj,nm)+miNew)/normby;
        end
    end
    for nj=1:smx
        for nm=1:NumberComp
            if isnan(mappedX(nj,nm))
               mappedX(nj,nm)=0; 
            end
            if isinf(mappedX(nj,nm))
                mappedX(nj,nm)=0;
            end
        end
    end    

    clusters=round(NumberComp/2);  %Twice the dimensionality
    %Write Data
    fill1=strcat(dataSet,'-D-',num2str(NumberComp));
    fil1=strcat(pathhTOS,fill1,'-MET-',methd,'-CLU-',num2str(clusters),'.bin');
    fid = fopen(fil1, 'w');
    for tmpy=1:smx
        fwrite(fid, mappedX(tmpy,:), 'float');
    end
    fclose(fid);
    %FIND REFERENCE POINTS
    %Clustering Section
    try
        [IDX,C]=kmeans(mappedX,clusters);
    catch ER
        clusters=round(clusters/2); 
        try 
            [IDX,C]=kmeans(mappedX,clusters);
        catch ER2
            clusters=round(clusters/2); 
            try 
                [IDX,C]=kmeans(mappedX,clusters);
            catch ER3
                clusters=round(clusters/2); 
                [IDX,C]=kmeans(mappedX,clusters);
            end
        end
    end
    %IDX contains the cluster the points belong
    %C contains the cluster locations
    %Sequential mode of reference points
    fil2=strcat(pathhTOS,dataSet,'-D-',num2str(NumberComp),'-MET-',methd,'-CLU-',num2str(clusters),'.ref');
    fid = fopen(fil2, 'w');
    for tmpy=1:clusters
        fprintf(fid,'%f ', C(tmpy,:));
        fprintf(fid,'\n');
    end
    fclose(fid);     
    
    
    smt=size(t_points,1);
    miNew=abs(min(min(t_points)));
    maNew=max(max(t_points));
    normby=maNew+miNew;
    for nj=1:smt
        for nm=1:NumberComp
            t_points(nj,nm)=(t_points(nj,nm)+miNew)/normby;
        end
    end
    for nj=1:smt
        for nm=1:NumberComp
            if isnan(t_points(nj,nm))
               t_points(nj,nm)=0; 
            end
            if isinf(t_points(nj,nm))
                t_points(nj,nm)=0;
            end
        end
    end        

    %Write Data
    fill1=strcat(dataSet,'-D-',num2str(NumberComp));
    fil1=strcat(pathhTOS,fill1,'-MET-',methd,'-CLU-',num2str(clusters),'.quer');
    fid = fopen(fil1, 'w');
    for tmpy=1:smt
        fwrite(fid, t_points(tmpy,:), 'float');
    end
    fclose(fid);    
    
    no_dimsEV=tst;
    skip=0;
    end
 end %Reduction Methods loop
 end % Dimensions Loop
end %Dim Reduction execute or not?
%%Actual section of dimensionality reduction
if DIM_RED_T==1
 for lp_n_d=1:n_DM  %Loop through the targeted dimensions
 no_dimsEV=dim_to(lp_n_d);
 % Loop for 6 dim red techniques
 for dmr_t=1:8
    skip=0;
    if dmr_t==1
        %PCA STUFF
        methd='PCA';
        clear mappedX;
        clear t_points;
        [mappedX, mappingPCA] = compute_mapping(tmpH, 'PCA', no_dimsEV);
        %Generate New test set mapping
        %t_points = out_of_sample(tmpH2, mappingPCA);
    end
    if dmr_t==2
        methd='SVD';
        clear U;
        clear S;
        clear V;
        [U,S,V] = svd(zscore(tmpH),'econ');
        variances = diag(S).^2 / (size(tmpH,1)-1);
        varExplained = 100 * variances./sum(variances);
        Variance=cumsum(variances)./sum(variances);    %base 100          
        clear mappedX;
        clear t_points;
        %Training        
        mappedX = tmpH * V(:,1:no_dimsEV);    
        %Test Data
        %t_points= tmpH2 * V(:,1:no_dimsEV);    
    end    
    if dmr_t==3
        methd='KernelPCA';
        clear mappedX;
        clear t_points;
        [mappedX, mappingKPCA] = compute_mapping(tmpH, methd, no_dimsEV);
        %t_points = out_of_sample(tmpH2, mappingKPCA);
    end
    if dmr_t==4
        %FactorAnalysis Stuff
        methd='FactorAnalysis';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingFA] = compute_mapping(tmpH, methd, no_dimsEV);
        %Generate New test set mapping
        %t_points = out_of_sample(tmpH2, mappingFA);
    end
    if dmr_t==5
        methd='LLE';
        clear mappedX;
        clear t_points;    
        try
            [mappedX, mappingLLE] = compute_mapping(tmpH, methd, no_dimsEV);
        catch
            msg=char(strcat('LEE-',num2str(no_dimsEV),' failed!'))
            skip=1;
        end          
        if skip ~=1
            %Generate New test set mapping
            mappedX = out_of_sample(tmpH, mappingLLE);
            %t_points = out_of_sample(tmpH2, mappingLLE);
        end
    end
    if dmr_t==6
        methd='Laplacian';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingLAP] = compute_mapping(tmpH, methd, no_dimsEV);
        mappedX = out_of_sample(tmpH, mappingLAP);
        %t_points = out_of_sample(tmpH2, mappingLAP);
    end
    if dmr_t==7
        methd='Isomap';
        clear mappedX;
        clear t_points;        
        [mappedX, mappingISO] = compute_mapping(tmpH, methd, no_dimsEV);
        %Generate New test set mapping
        mappedX = out_of_sample(tmpH1, mappingISO);
        %t_points = out_of_sample(tmpH2, mappingISO);
    end
    if dmr_t==8
        methd='LPP';
        clear mappedX;
        clear t_points;        
        try
            [mappedX, mappingLPP] = compute_mapping(tmpH, methd, no_dimsEV);
        catch
            msg=char(strcat('LPP-',num2str(no_dimsEV),' failed because number of dimensions exceeds number of samples'))
            skip=1;
        end
        %Generate New test set mapping
        %if skip ~=1
        %    t_points = out_of_sample(tmpH2, mappingLPP);
        %end
    end
    if skip ~= 1
    [n1,n2]= size(mappedX);
    tst=no_dimsEV;
    if no_dimsEV ~= n2    %Very few times the mapping returns less than the desired dimensions, this will fix this problem
        tst=no_dimsEV;    %Random times when something goes wrong
        no_dimsEV = n2;
    end
    NumberComp=no_dimsEV;

    %WRITE TO INDEX FILES
    %Take out the bad apples  (LLE has many!)
    smx=size(mappedX,1);
    miNew=abs(min(min(mappedX)));
    maNew=max(max(mappedX));
    normby=maNew+miNew;
    for nj=1:smx
        for nm=1:NumberComp
            mappedX(nj,nm)=(mappedX(nj,nm)+miNew)/normby;
        end
    end
    for nj=1:smx
        for nm=1:NumberComp
            if isnan(mappedX(nj,nm))
               mappedX(nj,nm)=0; 
            end
            if isinf(mappedX(nj,nm))
                mappedX(nj,nm)=0;
            end
        end
    end    

    clusters=round(NumberComp/2);  %Twice the dimensionality
    %Write Data
    fill1=strcat(dataSet,'-D-',num2str(NumberComp));
    fil1=strcat(pathhTOS,fill1,'-MET-',methd,'-CLU-',num2str(clusters),'THESIS.bin');
    fid = fopen(fil1, 'w');
    for tmpy=1:smx
        fwrite(fid, mappedX(tmpy,:), 'float');
    end
    fclose(fid);
    %FIND REFERENCE POINTS
    %Clustering Section
    try
        [IDX,C]=kmeans(mappedX,clusters);
    catch ER
        clusters=round(clusters/2); 
        try 
            [IDX,C]=kmeans(mappedX,clusters);
        catch ER2
            clusters=round(clusters/2); 
            try 
                [IDX,C]=kmeans(mappedX,clusters);
            catch ER3
                clusters=round(clusters/2); 
                [IDX,C]=kmeans(mappedX,clusters);
            end
        end
    end
    %IDX contains the cluster the points belong
    %C contains the cluster locations
    %Sequential mode of reference points
    fil2=strcat(pathhTOS,dataSet,'-D-',num2str(NumberComp),'-MET-',methd,'-CLU-',num2str(clusters),'THESIS.ref');
    fid = fopen(fil2, 'w');
    for tmpy=1:clusters
        fprintf(fid,'%f ', C(tmpy,:));
        fprintf(fid,'\n');
    end
    fclose(fid);     
    no_dimsEV=tst;
    skip=0;
    end
 end %Reduction Methods loop
 end % Dimensions Loop    
end
msg='Indexing Module Demo has been completed. Check your output folder for results'    
    