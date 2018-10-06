classdef ThermoClustering < handle
    
    properties (SetAccess = private)
        resolution,
        rawObj,
        msObj,
        msNames,
        defaultOcc,
        markers,
        referMZ,        
    end
    
    properties (Dependent)
        nFiles;
    end
    
    methods
        function obj = ThermoClustering(resolution,occ,m1,m2,referMz)
            obj.rawObj = {};
            obj.msObj = {};
            obj.msNames = {};
            
            if ~exist('resolution','var')
                obj.resolution = 0.01;
            else
                obj.resolution = resolution;
            end
            
            if ~exist('occ','var')
                obj.defaultOcc = 10;
            else
                obj.defaultOcc = occ;
            end
            if exist('m1','var') && exist('m2','var')
                obj.markers = {m1,m2};
            else
                obj.markers = {};
            end
            if exist('referMz','var')
                obj.referMZ = referMz;
            else
                obj.referMZ = [];
            end
        end
        
        function r = get.nFiles(obj)
            r = length(obj.msNames);
        end
        
        function addFiles(obj,names)
            preLen = obj.nFiles;
            L = length(names);          
            msFilePath = cell(L,1);
            for m = 1:L
                [fn,fp,index] = uigetfile('*.raw',sprintf('Select file for %s...',names{m}));
                if index
                    msFilePath{m} = strcat(fp,fn);
                else
                    fprintf(1,'failed to add %s\n',names{m});
                    msFilePath{m} = [];
                end
            end
            I = cellfun(@(x)isempty(x),msFilePath);
            names = names(I==0);
            L = sum(I==0);
            if L > 0
                obj.msNames((preLen+1):(preLen+L)) = names;
                obj.rawObj((preLen+1):(preLen+L)) = cell(L,1);
                obj.msObj((preLen+1):(preLen+L)) = cell(L,1);
                for m = 1:L
                    fprintf(1,'Loading %s ...\n',names{m});
                    obj.rawObj{preLen+m} = ThermoMSRaw(msFilePath{m},0);
                    obj.reParse(preLen+m);
                end
            end
        end
        
        function filterOcc(obj,n,varargin)
            if nargin == 1 || isempty(n)
                n = 1:obj.nFiles;
            end
            
            if nargin > 2
                cellfun(@(x)x.featureSelectByOcc(varargin{:}),obj.msObj(n),'UniformOutput',0);
            else               
                cellfun(@(x)x.featureSelectByOcc(obj.defaultOcc,0),obj.msObj(n),'UniformOutput',0);
            end
        end
        
        function filterAssem(obj,n,varargin)
            if isempty(obj.referMZ)
                error('reference MZ not set');
            end
            if nargin == 1 || isempty(n)
                n = 1:obj.nFiles;
            end
            
            if nargin > 2
                cellfun(@(x)x.featureSelectByAssemRatio(obj.referMZ,varargin{:}),...
                    obj.msObj(n),'UniformOutput',0);
            else
                cellfun(@(x)x.featureSelectByAssemRatio(obj.referMZ,0.05,0),...
                    obj.msObj(n),'UniformOutput',0);
            end
        end
        
        function doBLSubstrate(obj,n,varargin)
            if nargin == 1 || isempty(n)
                n = 1:obj.nFiles;
            end
            
            if nargin > 2
                cellfun(@(x)x.blSubstrate(varargin{:}),obj.msObj(n),'UniformOutput',0);
            else
                cellfun(@(x)x.blSubstrate(0.01,0.05),obj.msObj(n),'UniformOutput',0);
            end
        end
        
        function show(obj,wd)
            if ~exist('wd','var')
                wd = 100;
            end
            figure;
            for m = 1:obj.nFiles
                bar(subplot(obj.nFiles,1,m),obj.msObj{m}.mzList,mean(obj.msObj{m}.dataMat),wd);
                title(obj.msNames{m});
            end
        end
        
        function reParse(obj,n)
            if nargin == 1
                n = 1:obj.nFiles;
            end
            obj.msObj(n) = cellfun(@(x)AlignedMSSet(x,obj.resolution),obj.rawObj(n),'UniformOutput',0);
        end
        
        function cluster(obj,n,isCsv)
            if ~exist('isCsv','var')
                isCsv = 0;
            end
            if isempty(obj.markers)
                error('markers not set');
            end
            L = length(n);
            if ~all(and(n>=1,n<=obj.nFiles))
                error('exists invalid index');
            end
            res = cell(L,1);
            nCell = zeros(L,1);
            tag = [];
            for m = 1:L
                res{m} = obj.msObj{n(m)}.toSC(struct('r1',obj.markers{1},'r2',obj.markers{2},'csv',isCsv));
                res{m} = res{m}.mz;
                nCell(m) = size(res{m},1);
                tag = [tag;m*ones(nCell(m),1)];
            end
            mzLists = cellfun(@(x)x.mzList,obj.msObj(n),'UniformOutput',0);
            mz = mzLists{1};
            for m = 2:L
                mz = union(mz,mzLists{m});
            end
            mz = sort(mz);
            mat = zeros(sum(nCell),length(mzLists));
            for m = 1:obj.nFiles
                [~,I] = ismember(mzLists{m},mz);
                mat((1+sum(nCell(1:(m-1)))):sum(nCell(1:m)),I) = res{m};
            end
            
            xy = tsne(mat./max(mat,[],2));
            gscatter(xy(:,1),xy(:,2),tag);         
        end
        
        function setReferMZ(obj,r)
            obj.referMZ = r;
        end
        
        function setMarkers(obj,r1,r2)
            obj.markers = {r1,r2};
        end
        
        function r = nPeaks(obj)
            r = cellfun(@(x)length(x.mzList),obj.msObj);
        end
    end
    
end

