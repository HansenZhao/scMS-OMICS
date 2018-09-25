classdef SIMSImTxtRaw < handle
    
    properties (SetAccess=private)
        fileName, % interface
        scanLength,
        shotsPerPixel
        imageSize,
        instName,
        readMethod,
        mz,
        intens,
        additionAttr,
        scanID
    end
    
    properties (Dependent)
        fileNum
        scanNumber
    end
    
    methods
        function obj = SIMSImTxtRaw(folderName)
            if nargin == 0
                folderName = uigetdir();
            end
            
            obj.fileName = '';
            
            [obj.mz,obj.intens,obj.scanID] = deal({});                   
            obj.scanNumber = 0;
            obj.additionAttr = struct();
            obj.additionAttr.X = [];
            obj.additionAttr.Y = [];
            obj.addMSByFile(folderName);
        end            
       
        function res = get.fileNum(obj)
            if isempty(obj.fileName)
                res = 0;
            else
                res = 1 + length(strfind(obj.fileName,'|'));
            end
        end
        
        function res = get.scanNumber(obj)
            res = size(obj.intens,1);
        end
        
        function addMSByFile(obj,folderName)            
            [fns,fp] = listFile('*.txt',folderName);
            txtObj = SIMSTxtData(fp,fns{1});
            if obj.fileNum == 0
                obj.fileName = folderName;             
                obj.scanLength = txtObj.scanLength;
                obj.shotsPerPixel = txtObj.shotsPerPixel;
                obj.imageSize = txtObj.imageSize;
                obj.instName = 'SIMS-IMAGE';
            else
                if all(obj.scanLength == txtObj.scanLength &&...
                        obj.shotsPerPixel == txtObj.shotsPerPixel &&...
                        obj.imageSize == txtObj.imageSize)
                    obj.fileName = sprintf('%s | %s',obj.fileName, folderName);
                else
                    disp('Inconsistent file');
                    return;
                end               
            end
            
            sn = length(fns);
            offset = obj.scanNumber;
            charac = char(obj.fileNum+64);
            
            nMZ = length(obj.mz);
            obj.mz = [obj.mz,zeros(sn,1)];
            obj.intens = [obj.intens;zeros(obj.imageSize^2,length(obj.mz))];
            obj.scanID = [obj.scanID;cell(sn,1)];
            
            h = waitbar(0,'Parsing...');
            for m = 1:1:sn
                try
                    txtObj = SIMSTxtData(fp,fns{m});
                    if txtObj.mz > 0
                        I = ismember(txtObj.mz, obj.mz);
                        if I<0
                            obj.mz(nMZ+1) = txtObj.mz;
                            nMZ = nMZ + 1;
                        end
                        obj.intens(offset+(1:power(obj.imageSize,2)),I) = txtObj.rawMat(:);
                        obj.scanID{offset+m} = strcat(charac,num2str(m));
                    end
                catch
                    disp('a');
                end
                if mod(m,floor(sn/10)) == 0
                    waitbar(m/sn,h);
                end
            end
            [X,Y] = meshgrid(1:obj.imageSize);
            obj.additionAttr.X = [obj.additionAttr.X;X(:)];
            obj.additionAttr.Y = [obj.additionAttr.Y;Y(:)];
            
            obj.mz((nMZ+1):end) = [];
            obj.intens(:,(nMZ+1):end) = [];
            obj.scanID((nMZ+1):end) = [];
            [~,I] = sort(obj.mz);
            obj.mz = obj.mz(I);
            obj.intens = obj.intens(:,I);
            close(h);
        end
        
        function sn = sampleNumber(obj)
            % interface function
            sn = obj.scanNumber;
        end
        
        function [mzList,massIntens,attri] = getSample(obj,id)
            % interface function
            if isa(id,'double') && id > 0 && id <= obj.sampleNumber
                I = id;
            elseif isa(id,'char')
                [~,tmp] = ismember(id,obj.scanID);
                if tmp > 0
                    I = tmp;
                end
            else
                [mzList,massIntens,attri] = deal([]);
                return;
            end
            mzList = obj.mz;
            massIntens = obj.intens(I,:);
            attri = struct();
            attri.sampleID = obj.scanID{I};
            attri.x = obj.additionAttr.X(I);
            attri.y = obj.additionAttrY(I);
        end
        
        function res = getSourceID(obj)
            % interface function
            % 00-Thermo-RAW file
            % 01-SIMS-TXT file
            % isSpatial | isTimeSeries | hasParent
            res = '01100';
        end
    end
    
    methods (Static)
        function res = parsePyList(x)
            tmp = cell(x);
            if isempty(tmp)
                res = [];
                return;
            end
            if isa(tmp{1},'py.str')
                disp('s')
            else
                res = cellfun(@double,cell(x));
            end
        end
        
        function res = parsePyStr(x)
            res = char(x);
        end
    end
    
end

