classdef SparseMapArray < SecmImageArray
% SPARSEMAPARRAY Creates an array of sparse maps under clp-secm imaging system.
%   obj = SPARSEMAPARRAY(ticks,'random-location',radius,ndicts) Creates the map uniform
%   randomly with given dictionary radius and numbers.
% 
%   obj = SPARSEMAPARRAY(ticks,'random-all',radius,ndicts) Creates the map uniform
%   randomly with given dictionary radius and numbers and 
%   magnitude~N(1,1/2).
% 
%   obj = SPARSEMAP(ticks,'userdef',locations,magnitude) Create map by its 
%   input location and magnitude infos.
%
% SPARSEMAPARRAY methods:
%   SPY_MAP   - Show sparse pattern of map
%   
% SPARSEMAPARRAY public fields:
%   MAGNITUDE - magnitude of each dicts
%   LOCATIONS - (i,j,x,y) of dict locations
%   NDICTS    - number of dict
%
% For other methods, see SecmImage
%   
% See also SECMIMAGEARRAY
properties
    magnitude % vector(ndict,1); magnitude of each spike
    locations % vector(ndict,3); location of each spike
    ndicts    % scalar; number of dicts
end

methods
    function obj = SparseMapArray(ticks, k1, k2, varargin)
        % SPARSEMAP Construct the sparsity map
        gentype = varargin{1};
        if ~strcmp(gentype,'random-location') && ...
           ~strcmp(gentype,'random-all') && ...
           ~strcmp(gentype,'userdef')
            error('Wrong map generation type.');
        end

        % Assign input parameters
        obj   = obj@SecmImageArray(ticks,k1,k2);
        
        images = cell(k1,k2);
        
        for i = 1:k1
            for j = 1:k2
                images{i,j} = zeros(obj.nmeasures,obj.nmeasures);
            end
        end

        xylim = max(ticks);
        if strcmp(gentype,'random-location')
            radius = varargin{2};
            ndicts = varargin{3}; 
            magnitude = ones(ndicts,2);
        elseif strcmp(gentype,'random-all')
            radius = varargin{2};
            ndicts = varargin{3};
            magnitude = 1+randn(ndicts,1)/4;
        elseif strcmp(gentype,'userdef')
            locations = varargin{2};
            magnitude = varargin{3};
            ndicts = size(locations,1);
        end

        % Unif-random / Self assign the sparsity map
        idict = 1;
        while idict <= ndicts
            % Find location
            if strcmp(gentype,'random-location') || ...
               strcmp(gentype,'random-all') 
                r = (xylim-2*radius)*sqrt(rand());
                t = 2*pi*rand();
                [xloc,xidx] = closest(ticks,r*cos(t));
                [yloc,yidx] = closest(ticks,r*sin(t));
            elseif strcmp(gentype,'userdef')
                [xloc,xidx] = closest(ticks,locations(idict,1));
                [yloc,yidx] = closest(ticks,locations(idict,2));
            end
            
            i = randi(k1);
            j = randi(k2);
            
            % Update map
            if images{i,j}(yidx,xidx) == 0  % The location not assigned before
                images{i,j}(yidx,xidx) = magnitude(idict);
                locations(idict,:) = [i,j,xloc,yloc];
                idict = idict+1;
            elseif strcmp(gentype,'userdef') % The location is assigned
                error('Repeated locations');
            end
        end
        obj.images = images;
        obj.locations = locations;
        obj.ndicts = ndicts;
    end

    function downsample(obj, ~)
        % [TODO] Suggest implementing nearset neighbor interpolation.
        error('Sparsity map cannot be downsampled.');
    end

    % JW: I don't understand why this is here
    function sia = pos(obj)
        % Output and secmimage object
        sia = SecmImageArray(obj.ticks,obj.images);
        sia.pos();
    end

    function sia = soft(obj,lda)
        % Output and secmimage object
        sia = SecmImageArray(obj.ticks,obj.images);
        sia.soft(lda);
    end

    function spy_map(obj)
        % obj.SPY_MAP Draws sparsity pattern with coordinate system.
        gcf; 
        k = 0;
        
        for i = 1:size(obj.images,1)
            for j = 1:size(obj.images,2)
                k = k + 1;
                subplot(size(obj.images,1),size(obj.images,2),k);
                spy(obj.images{i,j});
                delta = floor(obj.nmeasures/10);
                idx = 1:delta:obj.nmeasures;
                set(gca,'Xtick',idx,...
                        'Xticklabel',obj.ticks(idx),...
                        'Ytick',idx,...
                        'Yticklabel',obj.ticks(idx),...
                        'YDir','normal');
                xlabel('Distance/mm')
                ylabel('Distance/mm')
            end
        end
    end
end
end

