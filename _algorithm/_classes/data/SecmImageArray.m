classdef SecmImageArray < SecmCoords
% SECMIMAGEARRAY Creates or modifies an array of images with SECM coordinate system. 
%   The images will always supported only on Ceff.
%
%   The basic object here is a cell array size k1 x k2 of images of the same size. An
%   SECMIMAGEARRAY interacts with another by matrix multiplication in a
%   kind of convolutional algebra: 
%
%    [A*B](i,j) = \sum_k A(i,k) conv B(k,j)
% 
%   obj = SECMIMAGEARRAY(ticks,k1,k2) Creates object with assigned coordinate
%   and empty image.
% 
%   obj = SECMIMAGEARRAY(ticks,images) Creates object with with assigned 
%   coordinates and input images. Here, images should be a cell array of
%   size k1 x k2
%
% SECMIMAGEARRAY is a handle object.
%   
% SECMIMAGEARRAY methods:
%   PLUS   -  Overload '+'; takes two SECMIMAGEARRAYs of the same size (k1, k2) and returns a new SECMIMAGEARRAY in which all of the images are summed.
%   MINUS  -  Overload '-'; same as plus, but with subtrction          
%   MTIMES -  Overload '*'; this acts on two SECMIMAGEARRAYs of compatible
%               size (k1,a), (a,k2) to produce a new SECMIMAGEARRAY whose i, j image is
%               given by the sum over a of the convolution of the (i,a) and (a,j)
%               images.
%   CTRANSPOSE - Overload '; Reverses each of the individual images and
%               conjugates the array. The implication here is that Z |-> D' * Z
%               is the adjoint of the map X |-> D * X
%   SOFT   -  Return obj with each of the images soft-thresholded by lda.
%   POS    -  Return obj with each of the images projected onto the
%               positive orthant
%   SUM    -  Return sum of image values.    
%   NORM   -  Return Frobenious norm of image.
%   MAX    -  Return max value in image.
%   INPROD -  Return inner product between image and input mask.
%   ROTATE -  Rotate image by an input angle counterclockwise.
%   SET_IMAGE    - Assign input image.
%   LINE_PROJECT - Project image to lines of given slope angles.
%   DOWNSAMPLE   - Decrease resolution of image.
%   DRAW_IMAGES   - Draw image with inherent coordinate system
%   
% SECMIMAGEARRAY public fields: 
%   IMAGES     - Discrete image with inherit coordinates. 
%   NMEASURES  - Number of pixels in x(y) directions.
%   TICKS      - The ticks in both xy direction in (mm)
%   RESOLUTION - Resolution of coordinates system in (mm)
%   XYLIM      - The boundary of coordinate in (mm) 
%
% See also SECMCOORDS, SCANLINES, DICTPROFILE, SPARSEMAP

properties
    images  % cell array of size k1 x k2 ... each entry is a matrix(nmeasures,nmeasures) representing an image.
end

methods
    function obj = SecmImageArray(ticks,varargin)
        % Constuct object SECMIMAGE of input coordinates.
        obj = obj@SecmCoords(ticks); 
        N = obj.nmeasures;
        switch nargin
            case 1; error('SECMIMAGEARRAY: not enough input arguments -- need to pass in either dimensions of the array, or a cell array of images');
            case 2; obj.images = varargin{1}; 
                k1 = size(obj.images,1);
                k2 = size(obj.images,2);
                
                for i = 1:k1
                    for j = 1:k2
                        obj.images{i,j} = (obj.Ceff).*obj.images{i,j};
                    end
                end
            case 3; 
                k1 = varargin{1}; k2 = varargin{2}; 
                obj.images = cell( k1, k2 );
                for i = 1:k1
                    for j = 1:k2
                        obj.images{i,j} = zeros(N,N);
                    end
                end
        end
    end
    


    function obj = plus(obj1,obj2)
    % obj = obj1 + obj2 Return a new SECMIMAGE with summed image 
        check_issame(obj1,obj2);
        images = cell(size(obj1.images));
        
        for i = 1:size(obj1.images,1)
            for j = 1:size(obj1.images,2)
                images{i,j} = obj1.images{i,j} + obj2.images{i,j};
            end
        end
        
        obj = SecmImageArray( obj1.ticks, images );    
    end
    
    function obj = minus(obj1,obj2)
    % obj = obj1 - obj2 Return a new SECMIMAGE with summed image 
        check_issame(obj1,obj2);
        images = cell(size(obj1.images));
        
        for i = 1:size(obj1.images,1)
            for j = 1:size(obj1.images,2)
                images{i,j} = obj1.images{i,j} - obj2.images{i,j};
            end
        end
        
        obj = SecmImageArray( obj1.ticks, images );    
    end
      
    function obj = times(obj1,obj2)
    % obj = obj1 .* obj2 Return a new SECMIMAGEARRAY with entrywise mult image.
        if isa(obj1,'SecmImageArray')
            obj = SecmImageArray(obj1.ticks,obj1.images);
            if isa(obj2,'SecmImageArray')
                check_issame(obj1,obj2);
                
                for i = 1:size(obj.images,1)
                    for j = 1:size(obj.images,2)
                        obj.images{i,j} = obj1.images{i,j} .* obj2.images{i,j};
                    end
                end
                
            else % obj2 is matrix
                
                for i = 1:size(obj.images,1)
                    for j = 1:size(obj.images,2)
                        obj.images{i,j} = obj1.images{i,j} .* obj2;
                    end
                end
            end
           
        elseif isa(obj1, 'cell') % obj1 is matrix array   %%%Shijie
            obj = SecmImageArray(obj2.ticks,obj2.images);
            
            for i = 1:size(obj.images,1)
                for j = 1:size(obj.images,2)
                    obj.images{i,j} = obj2.images{i,j} .* obj1{i,j};
                end
            end
                
        else % obj1 is matrix               
            obj = SecmImageArray(obj2.ticks,obj2.images);
            
            for i = 1:size(obj.images,1)
                for j = 1:size(obj.images,2)
                    obj.images{i,j} = obj2.images{i,j} .* obj1;
                end
            end
        end   
    end
    
    function obj = mtimes(obj1,obj2)
    % obj = obj1 * obj2 Return a new SECMIMAGEARRAY with convolution of images.
        if isa(obj1,'SecmImageArray')
            if isa(obj2,'SecmImageArray')
                
                % check compatibility 
                if ~isequal(obj1.ticks,obj2.ticks)
                    error('SecmImageArray operation available only in same coordinates')
                end
                if ~isequal(size(obj1.images,2),size(obj2.images,1))
                    error('Inner dimension needs to be the same for multiplication')
                end
                
                obj = SecmImageArray(obj1.ticks,size(obj1.images,1),size(obj2.images,2));
                
                for i = 1:size(obj.images,1)
                    for j = 1:size(obj.images,2)
                        for k = 1:size(obj1.images,2)
                            fimage = SecmCoords.fft(obj1.images{i,k}).* ...
                            SecmCoords.fft(obj2.images{k,j});
                            obj.images{i,j} = obj.images{i,j} + real(SecmCoords.ifft(fimage));
                        end
                    end
                end
            else % obj2 is scalar
                obj = SecmImageArray(obj1.ticks,size(obj1.images,1),size(obj1.images,2));
                
                for i = 1:size(obj.images,1)
                    for j = 1:size(obj.images,2)
                        obj.images{i,j} = obj1.images{i,j} * obj2;
                    end
                end
            end
        else % obj1 is scalar
            obj = SecmImageArray(obj2.ticks,size(obj2.images,1),size(obj2.images,2));
            
            for i = 1:size(obj.images,1)
                for j = 1:size(obj.images,2)
                    obj.images{i,j} = obj2.images{i,j} * obj1;
                end
            end
        end        
    end
           
    function obj = soft(obj,lda)
        if isa(lda,'cell') %%% Shijie: lda is matrix array
            for i = 1:size(obj.images,1)
                for j = 1:size(obj.images,2)
                    obj.images{i,j} = soft(obj.images{i,j},lda{i,j}); 
                end
            end
        else % lda is matrix
            for i = 1:size(obj.images,1)
                for j = 1:size(obj.images,2)
                    obj.images{i,j} = soft(obj.images{i,j},lda); 
                end
            end
        end
    end

    function obj = pos(obj)
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                obj.images{i,j} = max(obj.images{i,j},0); 
            end
        end
    end

    
    function s = sum(obj) 
        s = zeros(size(obj.images));
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                s(i,j) = sum(sum(obj.images{i,j}));
            end
        end
    end

    function s = norm(obj)
        s = 0; 
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                s = s + (norm(obj.images{i,j},'fro'))^2; 
            end
        end
        s = sqrt(s);
    end
    % s = obj.NORM(); Return Frobenious norm of image array.

    function s = max(obj)
        s = zeros(size(obj.images));
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                s(i,j) = max(max(obj.images{i,j})); 
            end
        end
    end
    % s = obj.MAX(); Return max value of image.

    function s = inprod(obj1,obj2) 
    % s = obj.INPROD(mask); Return inner product of two image arrays
        check_issame(obj1,obj2);
        s = 0; 
        
        for i = 1:size(obj1.images,1)
            for j = 1:size(obj2.images,2)
                s = s + sum(sum(obj1.images{i,j} .* obj2.images{i,j}));
            end
        end 
    end 
    

    function set_images(obj,imgs)
        obj.images = cell(size(imgs));
        for i = 1:size(imgs,1)
            for j = 1:size(imgs,2)
                obj.images{i,j} = (obj.Ceff).*imgs{i,j}; 
            end
        end
    end
    % obj.SET_IMAGES(img); Set image as img

    function obj = rotate(obj,angle)
    % obj = obj.ROTATE(angles); Rotate image by angle counterclockwise
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)           
                obj.images{i,j} = obj.Ceff.*(obj.rotate_image(obj.images{i,j},angle));
            end
        end
    end
    
    function objnew = ctranspose(obj)
        
        imgs = cell( size(obj.images,2), size(obj.images,1) );
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                imgs{j,i} = conj(fliplr(flipud(obj.images{i,j})));
            end
        end
        
        objnew = SecmImageArray(obj.ticks,imgs);
    end
    
    function downsample(obj,downrate)
    % obj.DOWNSAMPLE(downrate); Downsample image by integer rate. 
        idxlist = 1:downrate:length(obj.ticks);
        downsample@SecmCoords(obj,downrate); 
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                obj.images{i,j} = obj.Ceff.*(obj.images{i,j}(idxlist,idxlist));
            end
        end
    end  
    
    function lines = line_project(obj,params)
    % lines = obj.LINE_PROJECT(params) Line scan of SECM image.
    % lines = obj.LINE_PROJECT() Line scan of with clpconfig parameters. 
    %
    % The algorithm consists of four steps:
    % 1. Scan with CLP of slope (cos(t),sin(t)) of input angle t.
    %    After rotate the image by input angles, each scan starts from the
    %    largest y-index to the smallest.
    % 2. Shift each scanlines with input shifts.
    % 3. Entrywise multiply each scanlines with input intensity.
    % 4. Convolve each scanlines with input point spread function.
    %
    % LINE_PROJECT has an adjoint operator called back projection.
    %
    % See also SCANLINES.BACK_PROJECT, CLPCONFIG
        if nargin == 1
            params = ProbeParams('config',obj.ticks);
        end
        
        if ~isequal(size(obj.images),[1 1])
            error('Line Projection is only currently for SecmImageArray of size 1,1');
        end
        
        angles = params.angles.value;
        nlines = length(angles);
        currents = zeros(obj.nmeasures, nlines);

        % Line projection with assigned slope angles
        for I = 1:nlines
            currents(:,I) = flip(sum(obj.rotate_image(obj.images{1,1},angles(I)),2));                     
        end
        lines = ScanLines(obj.ticks, currents, params);

        % Shift each lines with assigned shifts
        if params.shifts.isactive
            s = params.shifts.func(params.shifts.value);
            lines.shift(s);
        end

        % Multiply the current data with specified intensity
        if params.intensity.isactive
            m = params.intensity.func(params.intensity.value);
            lines.mult(m);
        end

        % Convolve the point spread function
        if params.psf.isactive
            p = params.psf.func(params.psf.value);
            lines.conv(p);
            if norm(p) == 0
                error('Point spread function should be nonzero');
            end
        end  
    end
    
    function draw_image(obj)
        h = gca;
        x = obj.ticks;
        y = obj.ticks;
        imagesc(x,y,obj.images{1,1});
        set(h,'YDir','normal');
        xlabel('Distance/mm'); xtickformat('%.1f');
        ylabel('Distance/mm'); ytickformat('%.1f');
        c = colorbar(); c.Label.String = 'Currents/A';
                        c.FontSize = h.FontSize - 4;
                        c.Label.FontSize = h.FontSize;
    end
    
    function draw_images(obj)
    % obj.DRAW_IMAGES(); Draws image with coordinate system.
        gcf;
        k = 0;
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                k = k + 1;
                %subplot(size(obj.images,1),size(obj.images,2),k);   %%%Shijie
                figure(k); %%%Shijie
                h = gca;
                x = obj.ticks;
                y = obj.ticks;
                imagesc(x,y,obj.images{i,j});
                set(h,'YDir','normal');
                xlabel('Distance/mm'); xtickformat('%.1f');
                ylabel('Distance/mm'); ytickformat('%.1f');
                c = colorbar(); c.Label.String = 'Currents/A';
                                c.FontSize = h.FontSize - 4;
                                c.Label.FontSize = h.FontSize;
            end
        end
    end

    
    function draw_sum_images(obj)   %%% Shijie
    % obj.DRAW_SUM_IMAGES(); Draws sum image with coordinate system.
        gcf;
        k = 0;
        map_sum = zeros(size(obj.images{1,1}));
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                k = k + 1;
                h = gca;
                x = obj.ticks;
                y = obj.ticks;
                map_sum = map_sum + obj.images{i,j};  %%%Shijie: display the sum map
            end
        end
        imagesc(x,y,map_sum); 
        set(h,'YDir','normal');
        xlabel('Distance/mm'); xtickformat('%.1f');
        ylabel('Distance/mm'); ytickformat('%.1f');
        c = colorbar(); c.Label.String = 'Currents/A';
        c.FontSize = h.FontSize - 4;
        c.Label.FontSize = h.FontSize;
    end
    
    function imgarray = translate(obj, xdis, ydis)
    % obj.TRANSLATE(xdis,ydis)
    % xdis and ydis are the number of PIXELS to translate
    % xdis > 0 denotes translating to right-ward;
    % ydis > 0 denotes translating to up-ward.
        imgarray = obj;
        
        for i = 1:size(imgarray.images,1)
            for j = 1:size(imgarray.images,2)

                n = size(imgarray.images{i,j}, 1);
                if xdis > 0
                    imgarray.images{i,j} = [zeros(n,xdis), imgarray.images{i,j}(:,1:end-xdis)];
                elseif xdis < 0
                    imgarray.images{i,j} = [imgarray.images{i,j}(:,-xdis+1:end), zeros(n,-xdis)];
                end
                if ydis > 0
                    imgarray.images{i,j} = [imgarray.images{i,j}(ydis+1:end,:); zeros(ydis,n)];
                elseif ydis < 0
                    imgarray.images{i,j} = [zeros(-ydis,n); imgarray.images{i,j}(1:end+ydis,:)];
                end
            end
        end
    end
    
    function loss = square_loss(obj, img_truth)
    % obj.SQUARE_LOSS(img_truth); Calculates the square loss against ground
    % truth point-list.
    
        check_issame(obj,img_truth);
        
        loss = 0;
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)                
                loss = loss + sum(sum( (obj.images{i,j} - img_truth.images{i,j}).^2 ));
            end
        end
    end
    
    function loss = min_square_loss(obj, x_truth, D)
        
        error('JOHN: I have not implemented this for arrays yet.');
        % code below is copied from SecmImage
    % obj.MIN_SQUARE_LOSS(image_truth);
        TRANS_MAX = floor(size(obj.image,1) / 10);
        ROTAT_MAX = pi / 10;
        loss = inf;
        
        for i = -TRANS_MAX : TRANS_MAX
            for j = -TRANS_MAX : TRANS_MAX
                for r = linspace(-ROTAT_MAX, ROTAT_MAX, 10)
                    x = x_truth * [cos(r),-sin(r);sin(r),cos(r)];
                    x(:,1) = x(:,1) + i*obj.resolution;
                    x(:,2) = x(:,2) + j*obj.resolution;
                    
                    img = SecmImage(obj.ticks);
                    id = floor((x - obj.ticks(1)) / obj.resolution);
                    for k = 1:size(x,1)
                        img.image(id(k,1), id(k,2)) = 1;
                    end
                    img = img * D;

                    loss = min(loss, obj.square_loss(img));
                end
            end
        end
    end
    
    function normalize(obj)
    % obj.NORMALIZE(); Normalizes the images to have max value 1.
        m = max(max(obj.max()));
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                obj.images{i,j} = obj.images{i,j} / m;
            end
        end
    end
end

methods (Static)
    function img = rotate_image(img,angle)
    % img = obj.ROTATE_IMAGE(img,angle) Rotate image clockwise.
        img = fft_rotate(img,-angle);
    end
end

methods (Access = private)
    function check_issame(obj1,obj2)
        if ~isequal(obj1.ticks,obj2.ticks)
            error('SecmImageArray operation available only in same coordinates')
        end
        
        if ~isequal( size(obj1.images), size(obj2.images) )
            error('SecmImageArray operation requires arrays of the same size')
        end
    end
end

end % classdef