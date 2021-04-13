 % EXAMPLE5
%
%  Demonstrates reconstruction with a dictionary consisting of rotations of
%  a basic motif
close all;
initpkg

%% ===== Parameter Setting ====== %%
ticks = [-1:0.01:1];   % Distance tags of each line measurment (mm)
ndiscs = 8;           % Number of discs
disc_radius = 0.05;    % Disc radius (mm)
angles = [0:20:160]';  % Scan angles 


%% ===== Generate discs and ground truth image ====== %%
% Parameter of data, enable angles parameter only.
p        = ProbeParams(0);     % Initialize 0 parameters.
p.angles = ProbeParam(angles); % Set constant parameter angles.

% Generate equilateral triangle
L = 0.1; % triangle side length 
deltaTheta = (2*pi/3) / 40; % discretization in angle

trianglefunc = @(x,y) (y > -L/(2*sqrt(3))) & (sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))) & (-sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))); 
theta = 0:deltaTheta:(2*pi/3);
k = length(theta);
rotatedTriangleFunc = cell(1,k);

for i = 1:length(theta)
    rotatedTriangleFunc{i} = @(x,y) trianglefunc( x * cos(theta(i)) + y * sin(theta(i)), y * cos(theta(i)) - x * sin(theta(i)) );
end

D = DictProfileArray(ticks, rotatedTriangleFunc);

figure(1);
clf;
for i = 1:length(theta)
    imagesc(D.images{i});
    axis equal;
    pause(.1);
end

% Generate random map X0
X0 = SparseMapArray(ticks, k, 1, 'random-location', disc_radius, ndiscs);

% Generate simulated image Y
Y = D * X0;

% Generate lines R from image Y
R = Y.line_project(p);

% Generate back projection image LtR from lines R
LtR = R.back_project(); % -- really want to backproject to an SECM image

% Plot results
figure;
subplot(221); Y.draw_image();          title('True Image');
subplot(222); D.draw_image();          title('Kernel');
subplot(223); R.plot_lines([0,40,80]); title('Lines');
subplot(224); LtR.draw_image();        title('Back projection image');


%% ===== Reconstrct the image with IPalm algorithm package ====== %%
% Generate problem with formulation 'Calibrated Lasso'
lda = 0.5*max(max(max(pos(D'*D))));

prb = CalibLassoMultiMotif(R,D,p,lda);   %%%%


%pause;

% Generate algorithm 'IPalm' and solve the simulated problem.
alg = IPalmSecmSimul(prb,X0,D);
alg.set_maxiter(500);
figure; 
alg.solve(); 
 
