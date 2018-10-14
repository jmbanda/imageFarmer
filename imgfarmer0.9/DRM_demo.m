%% Dimensionality Reduction Module (DRM) Demo
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
% MORE INFO
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
index=segments*segments*parameters; %Number of elements in the Feature Vector
dataSet='Solar';                    %Dataset Being Manipulated
pathh='C:\Demo\';                   %Root folder of datasets
%Save too path
pathhTOS='C:\Demo\';                   %Root folder of datasets
exten='tif';
clasesNames={'Active Region' 'Coronal Jet' 'Emerging Flux' 'Filament' 'Filament Activation' 'Filament Eruption' 'Flare' 'Oscillation'};
exlen=length(char(exten))+1;
plt=0;       %0 for plots , 1 for no plots
%If there is need to re-run a section  (DEVELOPER MODE)  Default all should be 1
weka_write=1;                       %0 for no weka writting, 1 for weka writing
%Image Parameter List 
imgParam={'Entropy' 'Mean' 'StandardDev' 'FracDim' 'Moment3' 'Moment4' 'Uniformity' 'RelSm' 'TamDir' 'TamCont'};
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
%% Dimensionality Estimation (via PCA and SVD Components)
%This will produce an 8 element vector with the four PCA dimensional
%targets (from 96 to 99% of variance) and the four SVD dimensional targets
%(from 96 to 99% of variance)
%Example: target_dimensions=[ PCA # of components with 96% of variance,
%PCA # of components with 97% of variance, etc, etc]

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
    
%General Estimation Variables
target_dimensions=zeros(8,1); %this stores 8 target dimensions (4 for PCA, 4 for SVD)
dim_counter=1;

%Variance Calculation for PCA
clear COEFF;
clear SCORE;
clear latent;
clear tsquare;
[COEFF,SCORE,latent,tsquare] = princomp(tmpH);
Variance=cumsum(latent)./sum(latent);   
for thres=96:99
    for sm=1:(cells*parameters)
        if Variance(sm) >= (thres/100)
            target_dimensions(dim_counter)=sm;
            dim_counter=dim_counter+1;
            break
        end
    end    
end

%Plot the Variances for PCA
if plt==0   
    figure(2);
    set(2, 'Visible', 'off');
    plot(Variance,'DisplayName','Variance','YDataSource','Variance');
    tmp2='PCA Component Variance for ';
    name=char([strcat(pathhTOS,tmp2,'ALL ',' Dataset-',dataSet,'-',int2str(segments),'x',int2str(segments))]);
    tmpTLT2=char('PCA Component Variance for ');
    tmpTLT=char([strcat(tmpTLT2,'ALL',' Dataset-',dataSet)]);
    title(char(tmpTLT));        
    saveas(2,name,'jpg');
    close(2);    
    figure(3);
    set(3, 'Visible', 'off');
    variances=latent;
    %You can easily calculate the percent of the total variability
    %explained by each principal component.
    percent_explained = 100*variances/sum(variances);
    pareto(percent_explained)
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    tmp2='PCA percent of the total variability explained by each principal component for ';
    name=char([strcat(pathhTOS,tmp2,'ALL',' Dataset-',dataSet,'-',int2str(segments),'x',int2str(segments))]);
    tmpTLT2=char('PCA percent of the total variability explained by each principal component for ');
    tmpTLT=char([strcat(tmpTLT2,'ALL',' Dataset-',dataSet)]);
    title(char(tmpTLT));        
    saveas(3,name,'jpg');
    close(3);        
end

%Variance Calculation for SVD
clear U;
clear S;
clear V;
% V contains scores
[U,S,V] = svd(zscore(tmpH),'econ');
variances = diag(S).^2 / (size(tmpH,1)-1);
varExplained = 100 * variances./sum(variances);
Variance=cumsum(variances)./sum(variances);    %base 100    
for thres=96:99
    for sm=1:(cells*parameters)
        if Variance(sm) >= (thres/100)
            target_dimensions(dim_counter)=sm;
            dim_counter=dim_counter+1;
            break
        end
    end    
end

if plt==0
    figure(2);
    set(2, 'Visible', 'off');
    plot(Variance,'DisplayName','Variance','YDataSource','Variance');
    tmp2='SVD Comp Variance for ';
    name=char([strcat(pathhTOS,tmp2,'ALL',' Dataset-',dataSet,'-',int2str(segments),'x',int2str(segments))]);
    tmpTLT2=char('SVD Comp Variance for ');
    tmpTLT=char([strcat(tmpTLT2,'ALL',' Dataset-',dataSet)]);
    title(char(tmpTLT));        
    saveas(2,name,'jpg');
    close(2);    
    figure(3);
    set(3, 'Visible', 'off');
    %variances=latent;
    %You can easily calculate the percent of the total variability
    %explained by each principal component.
    percent_explained = 100*variances/sum(variances);
    pareto(percent_explained)
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    tmp2='SVD percent of the total variability explained by each principal component for ';
    name=char([strcat(pathhTOS,tmp2,'ALL',' Dataset-',dataSet,'-',int2str(segments),'x',int2str(segments))]);
    tmpTLT2=char('SVD percent of the total variability explained by each principal component for ');
    tmpTLT=char([strcat(tmpTLT2,'ALL',' Dataset-',dataSet)]);
    title(char(tmpTLT));        
    saveas(3,name,'jpg');
    close(3);        
end
%End of the dimensionality estimation
msg='Targeted dimensions for experimentation: '
target_dimensions
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
n_DM=8; %Number of dimensions
for lp_n_d=1:n_DM  %Loop through the targeted dimensions
    
no_dimsEV=target_dimensions(lp_n_d);

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
        [U,S,V] = svd(zscore(tmpH),'econ');
        variances = diag(S).^2 / (size(tmpH,1)-1);
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
    %%OUTPUT TO WEKA
    if skip ~= 1
    [n1,n2]= size(mappedX);
    tst=no_dimsEV;
    if no_dimsEV ~= n2
        tst=no_dimsEV;    %Random times when something goes wrong
        no_dimsEV = n2;
    end
    NumberComp=no_dimsEV;
    if weka_write==1
    %WEKA Writting of PROJ
    clear tempData;
    clear featureNames;   

    for images=1:sH1_1
        for components=1:NumberComp
            tempData(images,components)=num2cell(mappedX(images,components));
            tempData(images,NumberComp+1)=labelsData1(images);
        end
    end
    clear tmc1;
    variable=1;
    for components=1:NumberComp
        tmp=num2str(variable);
        featureNames(components) = {tmp};
        variable=variable+1;
    end
    featureNames(NumberComp+1)={'label'};
    classindex=NumberComp+1; %Columns
        
    tmc1=strcat('67-33-',dataSet,'-',methd,' Components Training N-',num2str(NumberComp),'-Feature-ALL-',int2str(segments),'x',int2str(segments));
    try 
        tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
    catch
        msg=char(strcat('Weka object creation Error on:',tmc1));
    end
    tmc5=strcat(pathhTOS,tmc1,'.arff');

    try
        saveARFF(tmc5,tempCon);     
    catch
        msg=char(strcat('Weka file creation Error on:',tmc5));
    end       

    
    %WEKA Writting FOR TEST SET
    clear tempData;
    clear featureNames;   

    for images=1:sH2_1
        for components=1:NumberComp
            tempData(images,components)=num2cell(t_points(images,components));
            tempData(images,NumberComp+1)=labelsData2(images);
        end
    end
    clear tmc1;
    variable=1;
    for components=1:NumberComp
        tmp=num2str(variable);
        featureNames(components) = {tmp};
        variable=variable+1;
    end
    featureNames(NumberComp+1)={'label'};
    classindex=NumberComp+1; %Columns
    
    tmc1=strcat('67-33-',dataSet,'-',methd,' Components Test N-',num2str(NumberComp),'-Feature-ALL-',int2str(segments),'x',int2str(segments));
    try 
        tempCon = matlab2weka(tmc1,featureNames,tempData,classindex);
    catch
        msg=char(strcat('Weka object creation Error on:',tmc1));
    end
    tmc5=strcat(pathhTOS,tmc1,'.arff');

    try
        saveARFF(tmc5,tempCon);     
    catch
        msg=char(strcat('Weka file creation Error on:',tmc5));
    end        
    no_dimsEV=tst;
    skip=0;
    end %Weka write skipping
  end
    
end %Reduction Methods loop
end % Dimensions Loop
msg='Dimensionality Reduction Module Demo has been completed. Check your output folder for results'