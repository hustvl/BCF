%function code
    %compile kdtree, your system should have installed a c++ compiler.
    mex('./kdtree/kdnn.cc');
    mex('./kdtree/kdtree.cc');
    
    %add paths to matlab.
    addpath(genpath('.'));

    % Read the sample image in
    im = imread('shapessm.jpg')
    if(size(im, 3) > 1)
       im = rgb2gray(im); 
    end
    
    % Find edges using the Canny operator with hysteresis thresholds of 0.1
    % and 0.2 with smoothing parameter sigma set to 1.
    edgeim = edge(im,'canny', [0.1 0.2], 1);
    %edgeim = im;
    % figure(1), imshow(edgeim); truesize(1);
    edgeim = im2bw(edgeim, 0.1);
    % Link edge pixels together into lists of sequential edge points, one
    % list for each edge contour.  Discard contours less than 10 pixels long.
    [edgelist, labelededgeim] = edgelink(edgeim, 10);
    
    % Display the labeled edge image with random colours for each
    % distinct edge in figure 2
    % drawedgelist(edgelist, size(im), 1, 'rand', 2); axis off       
    %relink two edges with the curvature at conjucton points with curvature
    %lower that 0.6.
    % 7 is the scale size to estimate the curvature.
    % 5 is the effective radius.(if the distance of the end points of two 
    % edges is less than 5, they should be taken into consideration.)
    % you should adjust the parameters to fit your system.
    edgelist =  relinkedge_iter(edgelist, 10, 5, 0.6);
    % Display the labeled edge image with random colours for each
    % distinct edge in figure 3    
    % drawedgelist(edgelist, size(im), 1, 'rand', 3); axis off
    %introduce break points at an edge with the curvature higher that 0.35
    % 7 is the scale size to estimate the curvature.
    % you should adjust the parameters to fit your system.    
    edgelist = edge_separate_iter(edgelist, 7, 0.3);
    % Display the labeled edge image with random colours for each
    % distinct edge in figure 4    
    % drawedgelist(edgelist, size(im), 1, 'rand', 4); axis off  
    
    
    
