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
        eleCounter,
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
                threshold = 0.001;
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
        
        function res = toSC(obj,options)
            if strcmp(obj.sourceID(1:2),'00')
                if ~and(isfield(options,'r1'),isfield(options,'r2'))                    
                    error('r1 and r2 should be provided');
                end
                if ~isfield(options,'method')
                    options.method = 'max';
                end
                res = obj.thermoSC(options);
            end
            
        end
        
        function newMZ = featureSelectByOcc(obj,occThres)
            occurence = full(sum(obj.dataMat>0));
            if ~exist('occThres','var')
                histogram(occurence,'BinWidth',20);
                answer = inputdlg('threshold','',1);
                if isempty(answer) || isnan(str2double(answer{1}))
                    newMZ = []; 
                    return;
                else
                    occThres = str2double(answer{1});
                end
            end
            I = occurence > occThres;
            newMZ = obj.mzList(I);
            answer = questdlg(sprintf('Commit %d -> %d ?',length(obj.mzList),length(newMZ)));
            if strcmp(answer,'Yes')
                obj.mzList = newMZ;
                obj.dataMat = obj.dataMat(:,I);              
            end
        end
        
        function newMZ = featureSelectByCluster(obj,maxClusterNum)       
            vec = zeros(maxClusterNum,1);
            for m = 1:maxClusterNum
                [~,~,sumd] = kmeans(full(obj.dataMat)',m,'Distance','correlation');
                vec(m) = sum(sumd);
            end
            figure('Position',[0,0,1000,400]);
            plot(subplot(121),1:maxClusterNum,vec);
%             xy = tsne(full(obj.dataMat)','Distance','correlation');
%             scatter(subplot(122),xy(:,1),xy(:,2),15,'filled');
            answer = inputdlg('cluster number','',1);
            if isempty(answer) || isnan(str2double(answer{1}))
                newMZ = [];
                return;
            else
                clusterNum = str2double(answer{1});
            end
            reptime = 10;
            summer = inf;
            c = 0;
            while(c<reptime)
                [tI,tC,tS] = kmeans(full(obj.dataMat)',clusterNum,'Distance','correlation');
                if sum(tS)<summer
                    summer = sum(tS);
                    I = tI;
                    C = tC;
                end
                c = c + 1;
            end
            t = cell2mat(obj.getFieldByName('sampleTime'));
            figure('Position',[0,0,400,1000]);
            for m = 1:clusterNum
                plot(subplot(clusterNum,1,m),t,C(m,:));
            end
            answer = inputdlg('selected cluster','',1);
            if isempty(answer) || isnan(str2double(answer{1}))
                newMZ = [];
                return;
            else
                clusterID = str2double(answer{1});
            end
            I2 = I==clusterID;
            answer = questdlg(sprintf('Commit %d -> %d ?',length(I2),sum(I2)));
            
            dist = pdist2(obj.dataMat(:,I2)',C(clusterID,:),'correlation');
            
            if strcmp(answer,'Yes')
                obj.mzList = obj.mzList(I2);
                obj.dataMat = obj.dataMat(:,I2);              
            end
            
            [~,I] = sort(dist);
            newMZ = obj.mzList(I);
                      
        end
    end
    
    methods (Access=private)
        function obj = parseThermoRaw(obj,dataHandle,accuarcy,threshold,combineMethod)
            
            obj.sourceFileName = dataHandle.fileName;
            massRange = dataHandle.massRange;
            obj.mzList = massRange(1):accuarcy:massRange(2);
            obj.capacity = 1000000;
            obj.eleCounter = 0;
            obj.dataMat = zeros(obj.capacity,3);
            
            [~,~,attri] = dataHandle.getSample(1);
            attriNames = fields(attri);
            L = length(attriNames);
            obj.attriMap = containers.Map();
            obj.attriCell = cell(L+1,1);
            for m = 1:L
                obj.attriMap(attriNames{m}) = m;
                obj.attriCell{m} = cell(dataHandle.sampleNumber(),1);
            end
            obj.attriMap('oldID') = L+1;
            obj.attriCell{L+1} = cell(dataHandle.sampleNumber(),1);
            obj.counter = 0;
            
            for m = 1:1:dataHandle.sampleNumber
                [mz,intens,attri] = dataHandle.getSample(m);
                
                if ~isempty(mz)
                    intens = intens/max(intens);
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
                        tmp = mz'==newMzList;
                        tmp2 = sum(tmp);
                        for h = 1:1:mzLength
                            if tmp2(h) > 1
                                newData(h) = func(intens(tmp(:,h)));
                            else
                                newData(h) = intens(tmp(:,h));
                            end
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
                    
                    if loc(1)<=0 && loc(1) >= -1
                        loc(1) = 1;
                    end
                                        
                    obj.addNewMS(loc,intens,attri);
                    obj.attriCell{obj.attriMap('oldID')}{obj.counter(1)} = dataHandle.scanID{m};
                end
                
                if mod(m,1000) == 0
                    fprintf(1,'%d/%d\n',m,dataHandle.sampleNumber);
                end
            end
            
            obj.dataMat((obj.eleCounter+1):end,:) = [];
            obj.dataMat = sparse(obj.dataMat(:,1),obj.dataMat(:,2),obj.dataMat(:,3));
            I = find(sum(obj.dataMat)==0);
            obj.dataMat(:,I) = [];
            obj.mzList(:,I) = [];
            
                                       
            for m = 1:L+1
                obj.attriCell{m}((obj.counter+1):end) = [];
            end
        end
        
        function addNewMS(obj,loc,intens,attri)
            obj.counter = obj.counter + 1;
            L = length(loc);
            
            if (L+obj.eleCounter) > obj.capacity
                [obj.dataMat,tmp] = AlignedMSSet2.extendCap(obj.dataMat,1);
                obj.capacity = tmp(1);
            end
            
            obj.dataMat((1:L)+obj.eleCounter,:) = [ones(L,1)*obj.counter,loc(:),intens(:)];
            obj.eleCounter = obj.eleCounter + L;
            
            if ~isempty(attri)   
                attriNames = fields(attri);
                L = length(attriNames);
                for m = 1:L
                    obj.attriCell{obj.attriMap(attriNames{m})}{obj.counter(1)} =  extractfield(attri,attriNames{m});
                end
            end
        end
        
        function res = thermoSC(obj,options)
            scanids = cytoLikeSelect(obj,options.r1,options.r2);
            t = obj.getFieldByName('sampleTime');
            SNcell = getSCSN(scanids);
            L = length(SNcell);
            fprintf(1,'Detect %d cells\n',length(SNcell));
            times = zeros(L,1);
            scanID = zeros(L,1);
            cellID = (1:L)';
            mz = zeros(L,length(obj.mzList));
            if strcmp(options.method,'max')
                func = @(x)max(x,[],1);
            else
                func = @(x)mean(x,1);
            end
            for m = 1:1:L
                sns = SNcell{m};
                times(m) = mean(t(sns));
                scanID(m) = mean(sns);
                mz(m,:) = func(obj.dataMat(sns,:));
                if std(mz(m,:)) == 0
                    disp('s');
                end
            end
            figure;
            plot(t,obj.getMZRange(options.r1)); hold on;
            scatter(times,mean(obj.getMZRange(options.r1))*ones(L,1),'filled');
            [fn,fp,index] = uiputfile('*.csv');
            if index
                HScsvwrite(strcat(fp,fn),[scanID,times,mz],cellID,...
                    strcat('cell id,scan id,time,',strrep(array2str(obj.mzList),'  ',',')));
            end
            res = table(cellID,scanID,times,mz);
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

