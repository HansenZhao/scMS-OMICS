classdef SIMSTxtData < handle
    
    properties(SetAccess = private)
        sampleName;
        mz;
        scanLength;
        shotsPerPixel;
        imageSize;
        rawMat;
        fp;
        fn;
    end
    
    methods
        function obj = SIMSTxtData(fp,fn)
            if nargin == 0
                [fn,fp] = uigetfile('*.txt');
            end          
            fpath = strcat(fp,fn);
            obj.fp = fp;
            obj.fn = fn;
            try
                fid = fopen(fpath);
                
                curStr = fgetl(fid);
                while ~isempty(curStr) && curStr(1)=='#'
                    if contains(curStr,'Source Filename:')
                        m = regexp(curStr,'Source Filename: (?<name>\S+)','names');
                        obj.sampleName = m.name;
                    elseif contains(curStr,'Source Interval:')
                        m = regexp(curStr,'Source Interval: (?<name>\S+) \S+','names');
                        obj.mz = str2double(m.name);
                        if isnan(obj.mz)
                            obj.mz = -1;
                        end
                    elseif contains(curStr,'Field of View')
                        m = regexp(curStr,'Field of View: (?<name>\S+) \S+','names');
                        obj.scanLength = str2double(m.name);
                    elseif contains(curStr,'Shots')
                        m = regexp(curStr,'Shots per pixel: (?<name>\S+)','names');
                        obj.shotsPerPixel = str2double(m.name);
                    elseif contains(curStr,'Image Size')
                        m = regexp(curStr,'Image Size: (?<name>\S+) \S+','names');
                        obj.imageSize = str2double(m.name);
                    end
                    curStr = fgetl(fid);
                end
                tmp = fscanf(fid,'%d %d %f');
                obj.rawMat = reshape(tmp(3*(1:(256*256))),[obj.imageSize,obj.imageSize]);
%                 obj.rawMat = zeros(obj.imageSize);
%                 curPoint = 1;
%                 while ~isnumeric(curStr)
%                     try
%                         %m = str2double(strsplit(curStr,' '));
%                         %obj.rawMat(curPoint) = m(3);
%                         
%                         %m = strsplit(curStr,' ');
%                         %obj.rawMat(curPoint) = str2double(m{3});
%                         
%                         %m = regexp(curStr,'\S+ \S+ (?<name>\S+)','names');
%                         %obj.rawMat(curPoint) = str2double(m.name);
%                         
%                         m = sscanf(curStr,'%d %d %f');
%                         obj.rawMat(curPoint) = m(end);
%                         
%                         curPoint = curPoint + 1;
%                     catch
%                     end
%                     curStr = fgetl(fid);
%                 end
%                 if curPoint ~= (obj.imageSize^2 + 1)
%                     warning('Parsing number not much!');
%                 end
                fclose(fid);
            catch e
                if exist('fid','var')
                    fclose(fid);
                end
                disp(e);
            end
        end
        function show(obj,varargin)
            imagesc(obj.rawMat); 
            if nargin == 1
                colormap('hot');
            else
                colormap(varargin{1});
            end
            title(sprintf('%s - m/z: %.2f',obj.sampleName,obj.mz),'Interpreter','none');
        end
    end
    
end

