classdef ThermoClustering < handle
    
    properties (SetAccess = private)
        resolution,
        rawObj,
        msObj,
        msNames,
        defaultOcc,
        markers,
        referMZ,
        toCellRes
    end
    
    properties (Dependent)
        nFiles;
    end
    
    methods
        function obj = ThermoClustering(resolution,occ,m1,m2,referMz)
            obj.rawObj = {};
            obj.msObj = {};
            obj.toCellRes = {};
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
                    fprintf(1,'set %s path %s\n',names{m},msFilePath{m});
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
                obj.toCellRes((preLen+1):(preLen+L)) = cell(L,1);
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
        
        function [clusterRes] = cluster(obj,n,isCsv)
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
            sampleNames = cell(L,1);
            for m = 1:L
                if isempty(obj.toCellRes{n(m)})
                    obj.toCellRes{n(m)} = obj.msObj{n(m)}.toSC(...
                        struct('r1',obj.markers{1},'r2',obj.markers{2},'csv',isCsv));
                end
                res{m} = obj.toCellRes{n(m)}.mz;
                nCell(m) = size(res{m},1);
                tag = [tag;m*ones(nCell(m),1)];
                sampleNames{m} = obj.msNames{n(m)};
            end
            mzLists = cellfun(@(x)x.mzList,obj.msObj(n),'UniformOutput',0);
            mz = mzLists{1};
            for m = 2:L
                mz = union(mz,mzLists{m});
            end
            mz = sort(mz);
            I = find((mz(2:end) - mz(1:(end-1)))<0.1*obj.resolution);
            mz(I+1) = [];
            mat = zeros(sum(nCell),length(mzLists));
            for m = 1:L
                [~,I] = ismember(mzLists{m},mz);
                loc = find(I==0);
                for h = 1:length(loc)
                    [~,I(loc(h))] = min(abs(mz-mzLists{m}(loc(h)))); 
                end
                mat((1+sum(nCell(1:(m-1)))):sum(nCell(1:m)),I) = res{m};
            end
            
            clusterRes = struct();
            clusterRes.mat = mat;
            clusterRes.nCell = nCell;
            clusterRes.mz = mz;
            clusterRes.tag = tag;
            clusterRes.sampleNames = sampleNames;
            
%             xy = tsne(mat./max(mat,[],2));
%             h = gscatter(xy(:,1),xy(:,2),tag); 
%             for m = 1:L
%                 h(m).DisplayName = obj.msNames{n(m)};
%             end
%             xticks([]); yticks([]); xlabel('t-SNE Dim 1'); ylabel('t-SNE Dim 2');
        end
        
        function delSCRes(obj,n)
            obj.toCellRes{n} = {};
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
    
    methods (Static)
        function res = calCellOccIndex(clusterRes)
            mat = clusterRes.mat;
            mz = clusterRes.mz;
            tag = clusterRes.tag;
            tagSet = unique(tag);
            nTag = length(tagSet);     
            nMZ = length(mz);
            res = zeros(nMZ,nTag+1);
            res(:,1) = mz(:);
            for m = 1:nMZ
                for n = 1:nTag
                    res(m,n+1) = sum(mat(tag==tagSet(n),m) > 0)/sum(tag==tagSet(n));
                end
            end           
        end
        function [hf,mzUsed] = procTSNE(clusterRes,range,varargin)
            % procTSNE(clusterRes,[]) equal to procTSNE(clusterRes,[],'all',-1)
            % procTSNE(clusterRes,[],'any',0.2)
            % procTSNE(clusterRes,[],'all',0.2)
            % procTSNE(clusterRes,[300,500],'all',0.2);
            if nargin == 2
                comd = 'all';
                thres = -1;
            else
                comd = varargin{1};
                thres = varargin{2};
            end
            
            if isempty(range)
                range = [-1,inf];
            end
            
            mat = clusterRes.mat;
            tag = clusterRes.tag;
            sNames = clusterRes.sampleNames;
            hf = figure('Position',[100,100,520,500]);
            L = length(unique(tag));
            mzOccList = ThermoClustering.calCellOccIndex(clusterRes);
            switch comd
                case 'all'
                    I = and(all(mzOccList(:,2:end)>thres,2),...
                        and(clusterRes.mz>range(1),clusterRes.mz<range(2))');
                case 'any'
                    I = and(any(mzOccList(:,2:end)>thres,2),...
                        and(clusterRes.mz>range(1),clusterRes.mz<range(2))');
            end
            fprintf(1,'Filter Out %d peaks, remain: %d\n',sum(~I),sum(I));
            if all(I==0)
                disp('no peaks left');
                return;
            end
            mzUsed = clusterRes.mz(I);
            xy = tsne(mat./max(mat,[],2));
            h = gscatter(xy(:,1),xy(:,2),tag);
            for m = 1:L
                h(m).DisplayName = sNames{m};
            end
            xticks([]); yticks([]); xlabel('t-SNE Dim 1'); ylabel('t-SNE Dim 2');
        end
    end
    
end

