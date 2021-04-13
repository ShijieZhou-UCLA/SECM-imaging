classdef DictProfileArray < SecmImageArray
% DICTPROFILEARRAY 
% Creates an array of dictionary profile under clp-secm imaging system.
%   obj = DICTPROFILEARRAY(ticks,func) Creates a dictionary element using input
%   funtion handle func under coordinate system depend on ticks.
%
%   Please check its method in object SecmImage
%
% See also SECMIMAGE   
properties 
    funcs
end

methods
    function obj = DictProfileArray(ticks, funcs)
        % Creat dictionaries with provided cell array of function handles.
        [k1, k2] = size(funcs);
        obj = obj@SecmImageArray(ticks, k1, k2);
        obj.funcs = funcs;
        obj.funcs_to_images();
    end
    function downsample(obj,downrate)
        % Downsample the image of dictionary 
        downsample@SecmCoords(obj,downrate);
        obj.funcs_to_image();
    end
    function secmImageArray = pos(obj)
        % Output and secmimage object
        secmImageArray = SecmImageArray(obj.ticks,obj.images);
        secmImageArray.pos();
    end
    function secmImageArray = soft(obj,lda)
        % Output and secmimage object
        secmImageArray = SecmImageArray(obj.ticks,obj.images);
        secmImageArray.soft(lda);
    end
end

methods (Access = private)
    function funcs_to_images(obj)
        % Draw the indicator function as a map using monte-carlo method 
        [k1,k2] = size(obj.funcs);
        obj.images = cell(k1,k2);
        N = obj.nmeasures;        
        for i = 1:k1
            for j = 1:k2
                
                image = zeros(N,N);
                ns = 50;
                for ix = 1:N
                    for iy = 1:N
                        x = obj.ticks(ix); 
                        y = obj.ticks(iy);
                        xs = obj.resolution * (rand(ns,1)-0.5) + x;
                        ys = obj.resolution * (rand(ns,1)-0.5) + y;
                        image(iy,ix) = sum(obj.funcs{i,j}(xs,ys))/ns;
                    end
                end
                
                obj.images{i,j} = image;
            end
        end
    end
end
end