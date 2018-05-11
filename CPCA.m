%% CPCA for fMRI - Matlab
% 

%% Get Data
% https://github.com/rzlim08/CPCA_Example
% 
% If you have git
% 
% git clone https://github.com/rzlim08/CPCA_Example
% 
% if not, you can download the zip file

%% Imports and Dependencies
% 
% First add the packages that we'll need. I'll try to use as few dependencies as possible here. 
% Later releases of matlab have built-in functions for this. 
% Jimmy Shen's toolbox also is a good option 
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% 
% 
% Dependencies
% Matlab 

% Change this path if you did not set up your download in a standard way
restoredefaultpath
addpath(genpath([pwd filesep 'dependencies2']));

%% Data
% We'll also need to the data, so make sure we have a reference to that too. 

%%
% We can use the Matlab 'dir' function to read in all the images. 
% The \*.img string reads in all the files that end in '.img'
% Change this path to point to your data
data_path = [pwd filesep 'example_data_Single_Subject'  filesep 'example_data_Single_Subject'];
image_path = [data_path filesep 's01'];
scans = dir([image_path filesep '*.img'])

%%
% We can see there are 214 scans for this subject. We now want to get a reference to all
% the images and store them in a way that matlab can read. 
scan = scans(1);
ss_path = [image_path filesep scan.name];
ss_hdr = load_nii(ss_path)
ss_img = ss_hdr.img;
figure;
imagesc(squeeze(ss_img(:,:,17)))

%%
% Ok, so we have one image loaded in, we'll need to load all the images. 
path_to_scans = {};
for i = 1:size(scans)
    path_to_scans{end+1} = [image_path filesep scans(i).name];
end
% create a cell array of a path to each scan
path_to_scans = path_to_scans';
% read in scan headers
scan_hdr = cellfun(@load_nii, path_to_scans, 'UniformOutput', 0);

%%
% So now we have a cell array storing the spm headers for all the images. 
% We want to read those in and append them together to form a 2 dimensional matrix.
% We can use Matlab's cellfun for this.
 read_and_reshape = @(im_hdr)(double(reshape(im_hdr.img, 1,[])));
 im_cell = cellfun(read_and_reshape, scan_hdr, 'UniformOutput', 0);
 brain_scans = cell2mat(im_cell);
size(brain_scans)

%% Masking 

%%
% There are infinite ways of creating an image mask. An mask is an image of zeros 
% and ones that denote the areas to be analyzed. 
% We will use a precreated mask for the purpose of this analysis, 
% but we can easily make our own masking function.
% For example, we can create a mask by thresholding the scans at the global mean. 
global_mean = mean(mean(brain_scans))
thresholded = brain_scans>global_mean;
% for each voxel, if every scan is above the 
% global mean, include the voxel in the analysis
mask_test = floor(sum(thresholded)/size(brain_scans,1)); 
m = reshape(mask_test, size(ss_hdr.img));
figure;
imagesc(m(:,:,17))

% find where all scans are one, and take those out of the brain scans
masked_im = brain_scans(:, find(mask_test));
size(masked_im)


%%
% For the purposes of this workshop, we'll use a precreated mask.
% read in and flatten mask
mask = read_and_reshape(load_nii([data_path filesep 'mask.img']));

%%
% To get our data matrix, we will select only the voxels that
% we've included in our analysis. We call our data matrix Z. 
Z = brain_scans(:, find(mask));
size(Z)

%%
% We will want to standardize this matrix. 
% PCA should be mean centered at the very least to return meaningful results. 
Z = zscore(Z, 1);
all(abs(mean(Z))<0.01)
all(abs(std(Z)-1)<0.01)

%%
% This gives us a standardized data matrix. 
% We will now create a design matrix to represent the timing information

%% Creating the Design Matrix (G Matrix)
% 
% Since we want to get the neural responses to the task, 
% we want to constain our analysis to the scans where we expect to see these responses. 
% We create our design matrix by inserting ones into the places where we expect to see signals. 
% This is a similar procedure to that used in SPM's first level analysis. 
Letters2=[55.269 79.074 89.123 118.569 123.266 138.011 163.450 179.190 189.245];
Letters4=[5.025 10.383 20.765 69.681 74.377 108.515 132.986 158.425 199.293];
Letters6=[15.740 25.789 40.540 143.369 148.727 168.808 173.832 183.887 194.269];
Letters8=[30.486 35.844 45.237 59.960 64.656 84.098 94.481 113.211 128.290];
onsets = {Letters2; Letters4; Letters6; Letters8};
onsets = cellfun(@ceil, onsets, 'UniformOutput', 0);
conditions = size(onsets,1)
bins = 8
G = zeros(size(Z,1), bins*conditions);
% For each condition
for j = 1:conditions
    cond_onsets = onsets{j};
    % For each onset
    for i = 1:size(cond_onsets,2)
        % Take the scan at the time of the onset, 
        % and the next 7 scans after, insert ones into these scans
        G(cond_onsets(i):cond_onsets(i)+bins-1, ...
            (j-1)*bins+1:(j-1)*bins+bins) = ...
            G(cond_onsets(i):cond_onsets(i)+bins-1, ...
            (j-1)*bins+1:(j-1)*bins+bins)|...
            eye(bins);
    end
end

%% won't work on non-linux machines, SPM not compiled for the platform
% a = calculate_hrf_shape([Letters2' ones(9, 1)*10], 214, 3);
% b = calculate_hrf_shape([Letters4' ones(9, 1)*10], 214, 3);
% c = calculate_hrf_shape([Letters6' ones(9, 1)*10], 214, 3);
% d = calculate_hrf_shape([Letters8' ones(9, 1)*10], 214, 3);
% plot([a b c])
% title('HRF shape estimation for design')
% xlabel('scans')
% ylabel('hrf')

%%
% Now we have a design matrix for each of the different conditions and onsets.
% Finally, we want to standardize the G matrix. 
G = zscore(G,1);

%% Regression
% 
% We will then regress the Z matrix onto the G matrix. 
% This will give use the betas for the regression. 

% This uses Matlab's mldivide function. 
% We can also use the normal equations, which seem to work faster in Matlab 2016a
C = G\Z;  

%% SVD/PCA

%%
% Now we can do the singular value decomposition (SVD). 
% SVD is a matrix decomposition technique that will yield equivalent results to PCA. 
% Normally, with multiple participants we'd have to append all the matrices from each 
% participant together, but here we only have a single participant. 
disp(['rank of Z = ', num2str(rank(Z))])
disp(['rank of GC = ' num2str(rank(G*C))]) % The rank of GC is equal to the rank of G. 
% Thus we can see this as projecting Z onto a lower dimensional space
[U D V] = svds(G*C, 3);
% Create predictor weights by regressing the eigenvectors onto the design matrix
% and multiplying by the sqrt of the number of scans
P = (G\U)*sqrt(size(G,1));

%%
% Now we can make very disappointing plots of the predictor weights, 
% since we only have one participant. 
p_conds = reshape(P(:, 1), [8, 4])
figure;
plot(p_conds)
xlabel('time bin')
ylabel('predictor weight')
legend('2 Letter','4 Letter', '6 Letter', '8 Letter')

p_conds = reshape(P(:, 2), [8, 4])
figure;
plot(p_conds)
xlabel('time bin')
ylabel('predictor weight')
legend('2 Letter','4 Letter', '6 Letter', '8 Letter')

p_conds = reshape(P(:, 3), [8, 4])
figure;
plot(p_conds)
xlabel('time bin')
ylabel('predictor weight')
legend('2 Letter','4 Letter', '6 Letter', '8 Letter')

%% Extras
% Exploring the Impact of Analysis Software on Task fMRI Results
% 
%  https://doi.org/10.1101/285585
