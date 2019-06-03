classdef ThermoMSRaw < handle
    
    properties (SetAccess=private)
        fileName, % interface
        scanRange,
        timeRange,
        massRange,
        instName,
        tolerence,
        readMethod,
        mz,
        intens,
        scanID,
        parentMS,
        scanTime,
        scanNumber,
        additionAttr
    end
    
    properties (Dependent)
        fileNum
    end
    
    methods
        function obj = ThermoMSRaw(fileName,readMethod)
            obj.fileName = '';
            obj.scanNumber = 0;
            
            [obj.mz,obj.intens,obj.scanID] = deal({});
            [obj.scanTime,obj.parentMS] = deal([]);
                      
            if readMethod == 0
                [resolution,PPM] = deal(cell({}));
                obj.additionAttr = struct();
                obj.additionAttr.resolution = resolution;
                obj.additionAttr.PPM = PPM;
            end
            
            obj.addMSByFile(fileName,readMethod);
            
        end
        
        function res = get.fileNum(obj)
            if isempty(obj.fileName)
                res = 0;
            else
                res = 1 + length(strfind(obj.fileName,'|'));
            end
        end
        
        function addMSByFile(obj,fileName,readMethod)
            % readMethod: 0 for Peaks; 1 for Centronized; 2 for Profile
            try
                pyObj = py.ThermoMSReader.read(fileName,readMethod);
            catch e                
                rethrow(e);
            end
            
            if obj.fileNum == 0
                obj.fileName = fileName;
                obj.scanRange = ThermoMSRaw.parsePyList(pyObj.prop{'scan_range'});
                obj.timeRange = ThermoMSRaw.parsePyList(pyObj.prop{'time_range'});
                obj.massRange = ThermoMSRaw.parsePyList(pyObj.prop{'mass_range'});
                obj.instName = ThermoMSRaw.parsePyStr(pyObj.prop{'inst_name'});
                obj.tolerence = ThermoMSRaw.parsePyStr(pyObj.prop{'tolerence'});
                obj.readMethod = ThermoMSRaw.parsePyStr(pyObj.prop{'read_method'});
            else
                if all(obj.massRange==ThermoMSRaw.parsePyList(pyObj.prop{'mass_range'})) &&...
                   strcmp(obj.instName,ThermoMSRaw.parsePyStr(pyObj.prop{'inst_name'})) &&...
                   strcmp(obj.readMethod,ThermoMSRaw.parsePyStr(pyObj.prop{'read_method'}))
                    obj.fileName = sprintf('%s | %s',obj.fileName, fileName);
                else
                    disp('Inconsistent file');
                    return;
                end               
            end
            
            sn = double(pyObj.prop{'scan_number'});
            offset = obj.scanNumber;
            charac = char(obj.fileNum+64);
            obj.scanNumber = obj.scanNumber + sn;
            
            obj.mz = [obj.mz;cell(sn,1)];
            obj.intens = [obj.intens;cell(sn,1)];
            obj.scanID = [obj.scanID;cell(sn,1)];
            obj.scanTime = [obj.scanTime;zeros(sn,1)];
            obj.parentMS = [obj.parentMS;zeros(sn,1)];
            
            h = waitbar(0,'Parsing...');
            for m = 1:1:sn
                try
                    obj.mz{offset+m} = ThermoMSRaw.parsePyList(pyObj.fields{'mz_list'}{m});
                    obj.intens{offset+m} = ThermoMSRaw.parsePyList(pyObj.fields{'intens_list'}{m});
                    obj.scanID{offset+m} = strcat(charac,num2str(double(pyObj.fields{'scan_number'}{m})));
                    obj.scanTime(offset+m) = double(pyObj.fields{'scan_time'}{m});
                    filterStr = char(pyObj.fields{'filter_info'}{m});
                catch
                    disp('a');
                end
                tmp = strfind(filterStr,'ms2');
                if isempty(tmp)
                    obj.parentMS(offset+m) = nan;
                else
                    obj.parentMS(offset+m) = str2double(filterStr((tmp+4):(strfind(filterStr,'@')-1)));
                end 
                if mod(m,floor(sn/10)) == 0
                    waitbar(m/sn,h);
                end
            end
            
            if readMethod == 0
                obj.additionAttr.resolution = [obj.additionAttr.resolution;cell(sn,1)];
                obj.additionAttr.PPM = [obj.additionAttr.PPM;cell(sn,1)];
                for m = 1:1:sn
                    obj.additionAttr.resolution{offset+m} = ThermoMSRaw.parsePyList(pyObj.fields{'resolution'}{m});
                    obj.additionAttr.PPM{offset+m} = ThermoMSRaw.parsePyList(pyObj.fields{'ppm'}{m});
                end
            end
            
            close(h);
        end
        
        function sn = sampleNumber(obj)
            % interface function
            sn = obj.scanNumber;
        end
        
        function [mzList,massIntens,attri] = getSample(obj,id)
            % interface function
            if isa(id,'double') && id > 0 && id <= obj.sampleNumber
                mzList = obj.mz{id};
                massIntens = obj.intens{id};
                attri = struct();
                attri.sampleID = obj.scanID{id};
                attri.sampleTime = obj.scanTime(id);
                attri.parentMS = obj.parentMS(id);
%                 if strcmp(obj.readMethod,'Peaks')
%                     attri.resolution = obj.additionAttr.resolution{m};
%                     attri.PPM = obj.additionAttr.PPM{m};
%                 end
            elseif isa(id,'char')
                [~,tmp] = ismember(id,obj.scanID);
                if tmp > 0
                    mzList = obj.mz{tmp};
                    massIntens = obj.intens{tmp};
                    attri = struct();
                    attri.sampleID = obj.scanID{tmp};
                    attri.sampleTime = obj.scanTime(tmp);
                    attri.parentMS = obj.parentMS(tmp);
                end
            else
                [mzList,massIntens,attri] = deal([]);
            end         
        end
        
        function res = getSourceID(obj)
            % interface function
            % 00-Thermo-RAW file
            % 01-SIMS-TXT file
            % isSpatial | isTimeSeries | hasParent
            if all(isnan(obj.parentMS))
                res = '00010';
            else
                res = '00011';
            end
        end
        
        function res = queryMSByLowResMS(obj,r)
            if length(r) == 1
                r = [r-1/power(10,getFracPlace(r)),r+1/power(10,getFracPlace(r))];
            end
            mzList = zeros(obj.sampleNumber,1);
            for m = 1:obj.sampleNumber
                I = and(obj.mz{m}>=r(1),obj.mz{m}<=r(2));
                if sum(I) > 0
                    tmpMZ = obj.mz{m}(I);
                    tmpIntens = obj.intens{m}(I);
                    [~,I] = max(tmpIntens);
                    mzList(m) = tmpMZ(I);
                else
                    mzList(m) = nan;
                end
            end
            mzList = mzList(~isnan(mzList));
            if isempty(mzList)
                fprintf(1,'Query mz: %.4f failed\n',mean(r));
                res = [];
                return;
            end
            [mzs,~,ia] = unique(mzList);
            [v,freq] = mode(ia);
            if freq == 1
                res = mean(mzList);
            else
                res = mzs(v);
            end
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

