function J = bwsmooth_2( I, params )
%BWSMOOTH is a smoothing method for binarized black and white (bw) images.
% Image Processing Toolbox is required for the operation.
%
% J = bwsmooth( I )
% J = bwsmooth( I, params )
%
% input arguments :
%
%   I - the input black and white (bw) image
%   params - parameters (details are described later)
%
% output arguments :
%
%   J - the output black and white (bw) image
%
% parameters :
%
%   - params.bz_thresh = 0 ([0,-2] is the adequate range)
%     is the threshold for the final binarization.
%
%   - params.as_scale = 1/4 ([1/4,1] is the adequate range)
%     is the scale reduction ratio of distance image.
%
%   - params.debug = {0,1}
%     if 1 is set, distance images are shown for debugging
%
%
%
% REFERENCES
% [1] K.Shirai, Y.Endo, A.Kitadai, S.Inoue, N.Kurushima, H.Baba, A.Watanabe, M.Nakagawa,
%     "Character shape restoration of binarized historical documents by smoothing via geodesic morphology,"
%     in Proc. of IAPR Intel. Conf. on Doc. Analysis and Recognit. (ICDAR), 2013.
%

%
% parameter check
%
params = parameter_set( params );

%
% Step 1. Distance transform (DT):
% Conversion from the binary image to the distance image.
%
D = bw2dist( I, params );

figure(1003), imshow( I );
title('After Edge Detection') 
figure(1004), imshow( 1./ (1 + exp(-0.1*[D]) ) );
title('After Distance Transform') 

%
% Step 2. Anisotropic smoothing (AS):
% Smoothing for the distance image.
%
D = smooth_dist( D, params );

% figure(1005), imshow( 1./ (1 + exp(-0.1*[D]) ) );
% title('After Smoothing') 

%
% Step 3. Binarization (BZ):
% Conversion from the distance image to the binary image.
%
J = dist2bw( D, params );

figure(1006), imshow( J );
title('After Binarization') 

end

%
%-------------------------------------------------------------------------
%

function params = parameter_set( params_ )

%
% default parameters
%
params.bz_thresh = 0;
params.as_scale = 1;

% please add other user parameters


%
% overwrite the default parameters
%
names = fieldnames( params_ );

for i = 1:length( names )
	name = names{ i };
	params = setfield( params, name, getfield( params_, name ) );
end

end


%
% -------------------------------------------------------------------------
%

function D = bw2dist( B, params )
%BW2DIST
%
% D = BW2DIST( B, params )
%
% The binarized image B is converted to the distance image D.
%

% character or background
B_posi = B > 0; % char
B_nega = ~B_posi; % bg

% distance from the boundary
D = bwdist( bwperim( B_posi ), 'quasi-euclidean' );
% D = bwdist( bwperim( B_posi ), 'euclidean' );

% distinguish char. and bg. by signs
D_posi =   B_posi .* D;
D_nega = -(B_nega .* D);

D = D_posi + D_nega;

end

%
% -------------------------------------------------------------------------
%

function D_ = smooth_dist( D, params )
%SMOOTH_DIST
%
% J = SMOOTH_DIST( I, params )
%
% This function smooths an input image I by using Tschumperle's anisotropic
% diffusion.
%

%
% parameters of Tschumperle's PDE (Anisotropic diffusion)
%
nite = 30; % number of iterations

% smoothing parameters for gradient vector calculation
wsize = 5; % window size of gaussian function
wsigma = 1.5; % standard deviation

% diffusion parameters
enhance = 20; % enchancement degree. (enhance > 0)
dt = 0.25; % time step

scale = params.as_scale;

D_ = imresize( D, scale, 'cubic' );
D_= impdevec( D_, wsize, wsigma, enhance, dt, nite );
D_ = imresize( D_, size(D), 'cubic' );


% for DEBUG -----------------------------
if params.debug
figure(1001), imshow( 1./ (1 + exp(-0.1*[D, D_]) ) );
title( 'Smoothing for a distance image. before (left) and after (right)' );

DD = [D, D_];
DD( DD>0 ) = 0;
DD = -DD;
DD = exp( -DD/ 5 );
figure(1002), imagesc( DD ); axis image;
title( 'before / after with color enhancement for displaying' );

end
% for DEBUG -----------------------------

end

%
% -------------------------------------------------------------------------
%

function B = dist2bw( D, params )
%DIST2BW
%
% B = DIST2BW( D, params )
%
% A distance image D is converted to the binary image B.
%

thresh = params.bz_thresh;

B = D > thresh;

end


%
% -------------------------------------------------------------------------
%
%
% -------------------------------------------------------------------------
%

function J = impdevec( I, wsize, wsigma, enhance, dt, nIter )
%
% edge preserving smoothing for the multi channel image
% using Coherence Enhancing Flow
%
% J = impdeimpdecef( I, wsize, wsigma, permeation, nIter )
%
% inputarguments
%
%	I : Multi channel image
%
%
% wsize, wsigma : smoothing parameters for gradient vector calculation
%
%   wsize : window size of gaussian function
%   wsigma : standard deviation
%
% enhance : enchancement degree. (enhance > 0)
% dt      : time step
%	nIter   : number of iterations
%
%
% output arguments
%
%	J : smoothed image
%
%
%
% REFERENCE
%
% D. Tschumperle, R. Deriche,
% "Vector-valued image regularization with PDE's: A Common framework for
% different applications," in Proc. of IEEE Intl. Conf. on
% Comput. Vision and Pattern Recognit. (CVPR), 2003.
%
%

[sy, sx, sc] = size( I );
J = I;

% derivation operator mask
Hx = 0.5 * [ -1, 0, 1 ];
Hy = Hx';

Hxx = [ 1 -2 1 ];
Hyy = Hxx';
Hxy = 0.25 * [ 1 0 -1; 0 0 0; -1 0 1 ];

% smoothing operator mask
K = fspecial( 'gaussian', wsize, wsigma );

for iItr = 1:nIter

	% gradient of an image
	[Ix, Iy] = gradient_( J, Hx, Hy );

	% tensor of gradients
	[G11, G12, G21, G22] = structure_tensor( Ix, Iy, K );

	% Calculate the direction and the magnitude of gradients
	[R_min, R_max, Vx, Vy] = gradv( G11, G12, G21, G22 );

	r_std = enhance / ( std( sqrt( R_min( : ) ) ) + eps(1) );
 	R_min = r_std * R_min;
 	R_max = r_std * R_max;

	% set parameters for the diffusion tensor
	R2 = 1 ./ ( 1 + R_min + R_max );
	R1 = sqrt( R2 );
	
	[D11, D12, D21, D22] = diffusion_tensor( R1, R2, Vx, Vy, sc );

	% calculate trace( DH )
	[ Ixx, Ixy, Iyx, Iyy ] = hessian_( J, Hxx, Hxy, Hyy );
	DH = ( D11 .* Ixx + D12 .* Ixy ) + ( D21 .* Iyx + D22 .* Iyy );

	J = J + dt * DH;

end

end

%
% -------------------------------------------------------------------------
%

function [Ix, Iy] = gradient_( I, Hx, Hy )
	Ix = imfilter( I, Hx, 'replicate' );
	Iy = imfilter( I, Hy, 'replicate' );
end

%
% -------------------------------------------------------------------------
%

function [ Ixx, Ixy, Iyx, Iyy ] = hessian_( I, Hxx, Hxy, Hyy )
	Ixx = imfilter( I, Hxx, 'replicate' );
	Ixy = imfilter( I, Hxy, 'replicate' );
	Iyx = Ixy;
	Iyy = imfilter( I, Hyy, 'replicate' );
end

%
% -------------------------------------------------------------------------
%

function [G11, G12, G21, G22] = structure_tensor( Ix, Iy, K )
	G11 = sum( Ix.*Ix, 3 );
	G12 = sum( Ix.*Iy, 3 );
	G22 = sum( Iy.*Iy, 3 );

	G11 = imfilter( G11, K, 'replicate' );
	G12 = imfilter( G12, K, 'replicate' );
	G21 = G12;
	G22 = imfilter( G22, K, 'replicate' );
end

%
% -------------------------------------------------------------------------
%

function [ D11, D12, D21, D22 ] = diffusion_tensor( R1, R2, Vx, Vy, sc )
	D11 = R2 .* Vx.*Vx + R1 .* Vy.*Vy;
	D12 = R2 .* Vx.*Vy - R1 .* Vy.*Vx;
	D22 = R2 .* Vy.*Vy + R1 .* Vx.*Vx; 

	D11 = repmat( D11, [ 1 1 sc ] );
	D12 = repmat( D12, [ 1 1 sc ] );
	D21 = D12;
	D22 = repmat( D22, [ 1 1 sc ] );
end

%
% -------------------------------------------------------------------------
%

%
% Calculate the direction and the magnitude of gradients
% by using eigen-value analysis
%
function [ R_min, R_max, Vx, Vy ]...
	= gradv( G11, G12, G21, G22 )

	a = G11 + G22;
	b = G11 - G22;
	c = G12 .* G21;
	d = sqrt( b.^2 + 4 * c );

	R_min = 0.5 * ( a - d );
	R_max = 0.5 * ( a + d );

	R_min = max( R_min, 0 );
	R_max = max( R_max, 0 );

	t_max = atan2( R_max - G11, G12 );
	Vx = cos( t_max );
	Vy = sin( t_max );
end







