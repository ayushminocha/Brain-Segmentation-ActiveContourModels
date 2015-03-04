
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SEGMENTATION
% % CS269 Final Project

% % Professor: Demetri Terzopoulos 

% % Edgar Rios Piedra 
% % Harsh Jain 
% % Ayush Minocha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % PROGRAM USAGE

% % Select one folder that contains the 3d image files to be analyzed.
% % The program wil run automatically and output a segmentation for the white matter and gray matter for every 
% % slice where there is enough tissue to make a segmentation

% % It is important that the Matalb that is running this program has the SPM and FSL libraries installed (added to path). 
% % Since we are making use of some of their utlities to preprocess the cerebral images

% % Reference for SPM: http://www.fil.ion.ucl.ac.uk/spm/
% % Reference for FSL: http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL

% % Thanks!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Brain_Segmentation()

n = 1;
Original_Path = pwd;

fprintf('Selecting folder where images are located \n');
path = uigetdir(); 

cd (path)
name = dir('*.img');
number_of_files = size(name,1);

for i=1:number_of_files;   
    names(1,i) = cellstr(name(i).name);  
end

save('Brain_Files.mat')

disp('--------------------------------------------------------------------------')
disp('Pre-processing brain images');
disp('--------------------------------------------------------------------------')
                                 
Preprocessing_SPM (names, path, number_of_files);                          
Preprocessing_FSL (names, number_of_files);                              
%save('TestAfterBET.mat')

disp('--------------------------------------------------------------------------')
disp('Processing brain images');
disp('--------------------------------------------------------------------------')

for i=1:number_of_files                                                  
    
    NameNoImg = strrep(names(i), '.img', '');
    InputImages_Name(1,i) = strcat('BET_', NameNoImg, '.nii.gz');    
    InputMasks_Name(1,i)  = strcat('BET_', NameNoImg, '_mask.nii.gz'); 
    
end

for i = 1:number_of_files
    
    fprintf('---- Creating initial mask for brain #%d of %d ---- \n',i,number_of_files);
    
    [InputImages_3D] = Get_Stack3D (InputImages_Name(1,i));
    [InputMasks_3D]  = Get_Stack3D (InputMasks_Name (1,i));
    [X,Y,Z] = size(InputImages_3D);

    %for w=1:Z; imshow(InputImages_3D(:,:,w),[]); pause; end               % Display
    
    [Initial_3DMaskWM, Initial_3DMaskGM] = Get3DMask(InputImages_3D, InputMasks_3D, X, Y, Z);   

    %WorkSpace = ['SnakesInput_' int2str(i)];
    %save(WorkSpace)

    fprintf('---- Refining initial mask for brain #%d of %d ---- \n',i,number_of_files);

    %Snakes part - Inputs --> Initial_3DMaskWM, Initial_3DMaskGM
    %              Can find this vairables on the mat files called SnakesInput_1.mat, SnakesInput_2.mat, ...
    %              There is one mat file for each of the 5 brains analyzed
    
    LesSnakes(Initial_3DMaskWM, Initial_3DMaskGM, BruinBrain);
    
end

disp('--------------------------------------------------------------------------')
disp('FINISHED EVERYTHING - THANKS FOR USING BRAINMAX2.0')
disp('--------------------------------------------------------------------------')

cd(Original_Path)

end 

function Preprocessing_SPM (names, path, number_of_files)  

%%%%% FILE HANDLING

for i=1:number_of_files
    path_img(1,i)    = strcat(path,'/',names(i),',1');
end

all_subjects (1,:) = path_img;                                    

path   = all_subjects(1,:);
path   = cellstr(path);

%%%%% MOTION CORRECTION

matlabbatch{1}.spm.spatial.realign.estwrite.data = {{path{1:number_of_files}}}';
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

%%%%% SEGMENTATION -SPM-

matlabbatch{2}.spm.spatial.preproc.channel.vols = {path{1:number_of_files}};
matlabbatch{2}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{2}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,1'};
matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,2'};
matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,3'};
matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,4'};
matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,5'};
matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm = {'~/MATLAB/spm12/tpm/TPM.nii,6'};
matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{2}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{2}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{2}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{2}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{2}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{2}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{2}.spm.spatial.preproc.warp.write = [0 0];

%%%%% REGISTRATION -SPM-

matlabbatch{3}.spm.spatial.coreg.estwrite.ref(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{3}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{3}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

spm('defaults', 'FMRI');                                                   
spm_jobman('serial', matlabbatch);
clear matlabbatch;

end

function Preprocessing_FSL (names, number_of_files)

disp('----------------')
disp('Running FSL BET');
disp('----------------')

% Prepare the file names

for i=1:number_of_files
    BET_input(1,i)  = strcat('rr',names(i));                              
    BET_output(1,i) = strcat('BET_',names(i));                           
end

for i=1:number_of_files
    SystemCommand(1,i) = strcat( {'setenv FREESURFER_HOME /Applications/freesurfer ; source $FREESURFER_HOME/SetUpFreeSurfer.csh ;  bet '} , BET_input(1,i) , {' '}, BET_output(1,i), ' -m -R' );   % -m to get brain mask and -R to run bet various times to get robust result
end

for i=1:number_of_files
    SystemTemp = SystemCommand{1,i};                                      
    system(SystemTemp);                                                    % SHELL call
end
clear SystemCommand   

end

function Images_3D = Get_Stack3D (Images_Name)  

Images_Name = Images_Name{1};               
Images_3D = load_untouch_nii(Images_Name);
Images_3D = Images_3D.img;
Images_3D = imrotate(Images_3D,90);
[X,Y,Z] = size(Images_3D);

% Load in sagital view, changing to axial view
Axial = zeros(Y,Z,X);
for i = 1:X
    Axial(:,:,i) = squeeze(Images_3D(i,:,:));
end
Images_3D = Axial; 
[X,Y,Z] = size(Images_3D);

Big = zeros(X*2,Y*2,Z*2);
for i=1:Z
    temp = Images_3D(:,:,i);
    temp = imresize(temp,2);
    Big(:,:,i) = temp;
end
Images_3D = Big;
[X,Y,Z] = size(Images_3D);

end

function [Initial_3DMaskWM, Initial_3DMaskGM ] = Get3DMask(InputImages_3D, InputMasks_3D, X, Y, Z)

    for i = 1:Z                                                           
        
        BruinBrain = InputImages_3D(:,:,i);

        % Analyzing slice only if at least 20% of the image has brain tissue (it is not just air) 
        Decision = size(nonzeros(BruinBrain))/(X*Y);                       
        Threshold = .2;                                                  

        if (Decision(1) > Threshold)
            
            % Histogram analysis to find possible gray matter and white matter regions
            Max = (max(max(BruinBrain)));
            BB = BruinBrain/Max;
            [Counts,Intensity] = imhist(BB);                               
            Intensity = Intensity * Max;
            Counts(1)=Counts(2);      

            for h = 4:(size(Intensity,1)-6)                            
                Initial = Counts(h-3);
                Final = Counts(h+3);
                Diff(h-3) = abs(Final-Initial);
            end

            [sDiff, index] = sort(Diff, 'descend');

            Idx1 = index(1);   Idx2 = index(2);
            Idx3 = index(3);   Idx4 = index(4);   Idx5 = index(5);

            Countidx1 = Diff(Idx1+3);   Countidx2 = Diff(Idx2+3);
            Countidx3 = Diff(Idx3+3);   Countidx4 = Diff(Idx4+3);   Countidx5 = Diff(Idx5+3);
            
            % avoiding zeros
            if (Countidx1 == 0); Countidx1=1; end
            if (Countidx2 == 0); Countidx2=1; end
            if (Countidx3 == 0); Countidx3=1; end
            if (Countidx4 == 0); Countidx4=1; end
            if (Countidx5 == 0); Countidx5=1; end
            
            % avoiding max
            if (Countidx1 > 256); Countidx1=256; end
            if (Countidx2 > 256); Countidx2=256; end
            if (Countidx3 > 256); Countidx3=256; end
            if (Countidx4 > 256); Countidx4=256; end
            if (Countidx5 > 256); Countidx5=256; end

            PossibleThresholds = [Intensity(Countidx1) Intensity(Countidx2) Intensity(Countidx3) Intensity(Countidx4) Intensity(Countidx5) ];
 
            [~, index2] = sort(PossibleThresholds, 'ascend');
            
            Pre_T1 = PossibleThresholds(index2(1));
            Pre_T2 = PossibleThresholds(index2(2));
            
            % Analyzing if cut values are not so close together
            
            Diffs = abs(Pre_T1-Pre_T2);
            if (Diffs <= 50 ); Pre_T2 = PossibleThresholds(index2(3)); end
            
            Diffs = abs(Pre_T1-Pre_T2);
            if (Diffs <= 50 ); Pre_T2 = PossibleThresholds(index2(4)); end
            
            Diffs = abs(Pre_T1-Pre_T2);
            if (Diffs <= 50 ); Pre_T2 = PossibleThresholds(index2(5)); end
            
            Threshold1 = Pre_T1;   % less --> GM  , more  --> WM
            Threshold2 = Pre_T2;   % less --> AIR , more --> GM

            Temp = InputMasks_3D;
            BruinBrain_Mask  = zeros(X,Y);
            BruinBrain_MaskW = zeros(X,Y);
            BruinBrain_MaskG = zeros(X,Y);
            
            for u = 1:X
                for v = 1:Y
                    pix = BruinBrain(u,v);                                                         
                    if(pix >= Threshold1 && pix <= (Threshold2-1)); BruinBrain_Mask(u,v) = 200; BruinBrain_MaskG(u,v) = 200; end    % gray matter
                    if(pix >= Threshold2); BruinBrain_Mask(u,v) = 400; BruinBrain_MaskW(u,v) = 400;end                              % white matter
                end
            end
            
            % Pruning initial results

            % Dilate, erode  and fill images
            
            SE = strel('disk',1); 
            BruinBrain_MaskW2 = imdilate(BruinBrain_MaskW,SE);
            BruinBrain_MaskW2 = imfill(BruinBrain_MaskW2,'holes');
            BruinBrain_MaskW2 = imerode(BruinBrain_MaskW2,SE);

            SE = strel('disk',2); 
            BruinBrain_MaskG2 = imdilate(BruinBrain_MaskG,SE);
            BruinBrain_MaskG2 = imfill(BruinBrain_MaskG2,'holes');
            SE = strel('disk',1);  
            BruinBrain_MaskG2 = imerode(BruinBrain_MaskG2,SE);

            %Remove small regions
            
            BruinBrain_MaskW3 = bwareaopen(BruinBrain_MaskW2, 50);
            BruinBrain_MaskW3 = BruinBrain_MaskW3 .* BruinBrain_MaskW2;
            %imshow(BruinBrain_MaskW3,[]);        
            BruinBrain_MaskG3 = bwareaopen(BruinBrain_MaskG2, 50);
            BruinBrain_MaskG3 = BruinBrain_MaskG3 .* BruinBrain_MaskG2;
            %imshow(BruinBrain_MaskG3,[]);
            
            MaskW = BruinBrain_MaskW3;
            MaskG = BruinBrain_MaskG3;

            Decision = size(nonzeros(MaskW))/(X*Y);                      
            Threshold = .10;      

            if (Decision(1) > Threshold)
                Initial_3DMaskWM (1,:,:,i)= MaskW;
            else
                Initial_3DMaskWM (1,:,:,i)= zeros(X,Y);
            end
            
            Decision = size(nonzeros(MaskG))/(X*Y);                      
            Threshold = .10; 
            
            if (Decision(1) > Threshold)
                Initial_3DMaskGM (1,:,:,i)= MaskG;
            else
                Initial_3DMaskGM (1,:,:,i)= zeros(X,Y);
            end
 
            fprintf('Finished initial segmentation for slice #%d of %d \n',i,Z);
            
        else
            
            % Code to run when skipping slice
            
            Initial_3DMaskWM (1,:,:,i)= zeros(X,Y);              
            Initial_3DMaskGM (1,:,:,i)= zeros(X,Y); 

            fprintf('Skipped slice #%d of %d, not enough brain tissue \n',i,Z);
            
        end

    end
        
end

function LesSnakes (Initial_3DMaskWM, Initial_3DMaskGM, BruinBrain)
getBoundary(Initial_3DMaskWM, Initial_3DMaskGM, BruinBrain);
end

