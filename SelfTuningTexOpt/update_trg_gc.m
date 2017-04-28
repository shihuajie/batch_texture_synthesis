function [ trg_gc ] = update_trg_gc(num_iter)
    addpath('CCA');
    %% read images
    dir = ['Results\starry\iter_' num2str(num_iter)];
    trg = imread([dir filesep 'img.png']);
    
    %% parameters setting and initial
    patchSize = 10;
    numLabel = 16;
    numBins = 16;
    data_trg = [];
    [rows cols channels] = size(trg);
 
    %% feature of texture, there are total 6 features can be use
    % 1. gabor filter bank to make a feature space
%     [ G_src, G_trg ] = GaborFilter( src, trg );
%     srcGabor_Color = gaborColorStat(G_src, src, patchSize, mask);
%     trgGabor_Color = gaborColorStat(G_trg, trg, patchSize, []);
%     disp('.......generator gabor filter bank OK !');
%     % add Gabor filter
%     data_src = cat(2, data_src, srcGabor_Color);
%     data_trg = cat(2, data_trg, trgGabor_Color);
%     disp('.......add gabor filter bank OK !');
    
    % 2. generator SalMap
    trgSal = SalMap(trg);
    trgSalDat = histfeat( trgSal, patchSize, numBins, [] );
    data_trg = cat(2, data_trg, trgSalDat);
    disp('.......add salMap OK !');

    % 5. generator the MR8 filter
    X = genMR8Textons(trg);
    t = tic;
    opts = statset('Display', 'final');
    [idx,ctrs] = kmeans(X,numLabel,'Options',opts);
    fprintf('Kmeans in %f seconds.\n', toc(t));
    labelIm = reshape(idx,rows,cols);
    trgTextonDat = histfeat( labelIm, patchSize, numBins, [] );
    data_trg = cat(2, data_trg, trgTextonDat);
    disp('add MR8 textons channel OK !');
    
    % 6. statistic the histogram feature of image
    trgImgHist = ImHistStat(trg, patchSize, numBins, []);
    data_trg = cat(2, data_trg, trgImgHist);
    disp('.......add the histogram of image OK !');

%% CCA
    load 'CCADat.mat';
    X = data_trg;
    X = X-repmat(mean(X),(cols-patchSize+1) * (rows-patchSize+1),1);
    trg_gc = mat2gray(reshape(X*A, rows-patchSize+1, cols-patchSize+1));
    trg_gc = imresize(trg_gc,size(trg(:,:,1)),'lanczos3');
    imwrite(trg_gc, 'Results\starry\trg_gc.png');
    imwrite(trg_gc, [dir '\trg_gc.png']);
    
end
