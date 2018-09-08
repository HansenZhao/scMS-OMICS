classdef AlignedMSSet < handle
   
    properties (SetAccess = private)
        sourceFileName,
        sourceID,
        accuracy,
        mzList,
        dataMat
    end
    
    properties (Access = private)
        capacity,
        attriMap,
        attriCell,
        counter    
    end
    
    properties (Dependent)
        sampleNumber
    end
    
    methods
        function obj = AlignedMSSet(h,accuarcy,threshold,combineMethod)
            if ~exist('accuarcy','var')
                accuarcy = 0.001;
            end
            if ~exist('threshold','var')
                threshold = 0.01;
            end
            if ~exist('combineMethod','var')
                combineMethod = 'max';
            end
            tmp = h.getSourceID();
            obj.sourceID = tmp;
            obj.accuracy = accuarcy;
            if strcmp(tmp(1:2),'00')
                obj = obj.parseThermoRaw(h,accuarcy,threshold,combineMethod);
            end
        end
        
        function [res,attris] = getMSByScanNumber(obj,sn)
            I = find(obj.attriCell{obj.attriMap('oldID')}==sn);
            if isempty(I)
                disp('Scan number not found');
            else
                res = obj.dataMat(I,:);
            end
            fn = obj.fields;
            L = length(fn);
            attris = struct();
            for m = 1:L
                attris = setfield(attris,fn{m},obj.attriCell{obj.attriMap(fn{m})}(I));
            end
        end
        
        function res = getMZRange(obj,r,method)
            if ~exist('method','var')
                method = 'max';
            end
            if length(r) == 1
                [v,I] = min(abs(obj.mzList-r));
                if v <= (obj.accuracy*0.5)
                    res = double(obj.dataMat(:,I));
                end
            elseif length(r) == 2
                I = and(obj.mzList>=r(1),obj.mzList<=r(2));
                res = obj.dataMat(:,I);
                if strcmp(method,'max')
                    res = max(res,[],2);
                else
                    res = mean(res,2);
                end
            end      
        end
        
        function res = fields(obj)
            res = obj.attriMap.keys;
        end
        
        function res = getFieldByName(obj,fn)
            if obj.attriMap.isKey(fn)
                res = obj.attriCell{obj.attriMap(fn)};
            end
        end
        
    end
    
    methods (Access=private)
        function obj = parseThermoRaw(obj,dataHandle,accuarcy,threshold,combineMethod)
            
            obj.sourceFileName = dataHandle.fileName;
            massRange = dataHandle.massRange;
            obj.mzList = massRange(1):accuarcy:massRange(2);
            obj.dataMat = sparse(dataHandle.scanNumber,length(obj.mzList));
            
            [~,~,attri] = dataHandle.getSample(1);
            attriNames = fields(attri);
            L = length(attriNames);
            obj.attriMap = containers.Map();
            obj.attriCell = cell(L+1,1);
            for m = 1:L
                obj.attriMap(attriNames{m}) = m;
                obj.attriCell{m} = zeros(dataHandle.sampleNumber(),1);
            end
            obj.attriMap('oldID') = L+1;
            obj.counter = 0;
            
            for m = 1:1:dataHandle.sampleNumber
%                 if m == 491
%                     disp('s');
%                 end
                [mz,intens,attri] = dataHandle.getSample(m);
                
                if ~isempty(mz)
                    if strcmp(dataHandle.readMethod,'Profile')
                        intens = interp1(mz,intens,obj.mzList);
                        mz = obj.mzList;
                    else
                        mz = round(mz/accuarcy)*accuarcy;
                        newMzList = unique(mz);
                        mzLength = length(newMzList);
                        newData = zeros(1,mzLength);
                        if strcmp(combineMethod,'mean')
                            func = @(x)mean(x);
                        else
                            func = @(x)max(x);
                        end
                        for h = 1:1:mzLength
                            newData(h) = func(intens(mz == newMzList(h)));
                        end
                        mz = newMzList; intens = newData;
                    end
                    
                    if threshold < 1
                        I = intens > (threshold * max(intens));
                    else
                        I = intens > threshold;
                    end
                    
                    intens = intens(I);  mz = mz(I);
                    loc = round((mz-massRange(1))/accuarcy)+1;
                                        
                    obj.addNewMS(loc,intens,attri);
                    obj.attriCell{obj.attriMap('oldID')}(obj.counter(1)) = m;
                end
                
                if mod(m,1000) == 0
                    fprintf(1,'%d/%d\n',m,dataHandle.sampleNumber);
                end
            end
            
            obj.dataMat((obj.counter+1):end,:) = [];
            I = find(sum(obj.dataMat)==0);
            obj.dataMat(:,I) = [];
            obj.mzList(:,I) = [];
            
                                       
            for m = 1:L+1
                obj.attriCell{m}((obj.counter+1):end) = [];
            end
        end
        
        function addNewMS(obj,loc,intens,attri)
%             [~,I] = ismember(mz,obj.mzList);
            obj.counter = obj.counter + 1;
%             tmp = find(I<=0);
%             if ~isempty(tmp)
%                 L = length(tmp);
%                 fprintf(1,'scan: %d, found: %.4f\n',n,mz(tmp(1)));             
%                 for m = 1:1:L
%                     I(tmp(m)) = find(abs(obj.mzList-mz(tmp(m)))<0.000001);
%                 end
%             end
            obj.dataMat(obj.counter,loc) = intens;
            
            if ~isempty(attri)   
                attriNames = fields(attri);
                L = length(attriNames);
                for m = 1:L
                    obj.attriCell{obj.attriMap(attriNames{m})}(obj.counter(1)) = extractfield(attri,attriNames{m});
                end
            end
        end
        
    end
    
    methods (Static)
        function intens = interpMASSIntens(mz,x,y)
            intens = interp1(x,y,mz,'linear');
        end
        function [newr,s] = extendCap(r,dim)
            tmp = size(r);
            if all(dim==1)
                if iscell(r)
                    newr = cell(tmp(1)*2,tmp(2));
                else
                    newr = zeros(tmp(1)*2,tmp(2));
                end
                newr(1:tmp(1),:) = r;
                s = size(newr);
            elseif all(dim==2)
                if iscell(r)
                    newr = cell(tmp(1),tmp(2)*2);
                else
                    newr = sparse(tmp(1),tmp(2)*2);
                end
                newr(:,1:tmp(2)) = r;
                s = size(newr);
            else
                if iscell(r)
                    newr = cell(tmp*2);
                else
                    newr = zeros(tmp*2);
                end          
                newr(1:tmp(1),1:tmp(2)) = r;
                s = size(newr);
            end
        end
    end
    
end

