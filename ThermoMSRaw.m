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
    
    methods
        function obj = ThermoMSRaw(fileName,readMethod)
            % readMethod: 0 for Peaks; 1 for Centronized; 2 for Profile
            try
                pyObj = py.ThermoMSReader.read(fileName,readMethod);
            catch e                
                rethrow(e);
            end
            obj.fileName = fileName;
            obj.scanRange = ThermoMSRaw.parsePyList(pyObj.prop{'scan_range'});
            obj.timeRange = ThermoMSRaw.parsePyList(pyObj.prop{'time_range'});
            obj.massRange = ThermoMSRaw.parsePyList(pyObj.prop{'mass_range'});
            obj.instName = ThermoMSRaw.parsePyStr(pyObj.prop{'inst_name'});
            obj.tolerence = ThermoMSRaw.parsePyStr(pyObj.prop{'tolerence'});
            obj.readMethod = ThermoMSRaw.parsePyStr(pyObj.prop{'read_method'});
            obj.scanNumber = double(pyObj.prop{'scan_number'});
            
            [obj.mz,obj.intens] = deal(cell(obj.scanNumber,1));
            [obj.scanID,obj.scanTime,obj.parentMS] = deal(zeros(obj.scanNumber,1));
            
            h = waitbar(0,'Parsing...');
            for m = 1:1:obj.scanNumber
                obj.mz{m} = ThermoMSRaw.parsePyList(pyObj.fields{'mz_list'}{m});
                obj.intens{m} = ThermoMSRaw.parsePyList(pyObj.fields{'intens_list'}{m});
                obj.scanID(m) = double(pyObj.fields{'scan_number'}{m});
                obj.scanTime(m) = double(pyObj.fields{'scan_time'}{m});
                filterStr = char(pyObj.fields{'filter_info'}{m});
                tmp = strfind(filterStr,'ms2');
                if isempty(tmp)
                    obj.parentMS(m) = nan;
                else
                    obj.parentMS(m) = str2double(filterStr((tmp+4):(strfind(filterStr,'@')-1)));
                end 
                if mod(m,floor(obj.scanNumber/10)) == 0
                    waitbar(m/obj.scanNumber,h);
                end
            end
            
            if readMethod == 0
                [resolution,PPM] = deal(cell(obj.scanNumber,1));
                for m = 1:1:obj.scanNumber
                    resolution{m} = ThermoMSRaw.parsePyList(pyObj.fields{'resolution'}{m});
                    PPM{m} = ThermoMSRaw.parsePyList(pyObj.fields{'ppm'}{m});
                end
                obj.additionAttr = struct();
                obj.additionAttr.resolution = resolution;
                obj.additionAttr.PPM = PPM;
            end
            
            close(h);
        end
        
        function sn = sampleNumber(obj)
            % interface function
            sn = obj.scanNumber;
        end
        
        function [mzList,massIntens,attri] = getSample(obj,id)
            % interface function
            if id > 0 && id <= obj.sampleNumber
                mzList = obj.mz{id};
                massIntens = obj.intens{id};
                attri = struct();
                attri.sampleID = obj.scanID(id);
                attri.sampleTime = obj.scanTime(id);
                attri.parentMS = obj.parentMS(id);
%                 if strcmp(obj.readMethod,'Peaks')
%                     attri.resolution = obj.additionAttr.resolution{m};
%                     attri.PPM = obj.additionAttr.PPM{m};
%                 end
            else
                [mzList,massIntens,attri] = deal([]);
            end         
        end
        
        function res = getSourceID(obj)
            % interface function
            % 00-Thermo-RAW file
            % 01-SIMS-TXT file
            % isSpatial | is
            if all(isnan(obj.parentMS))
                res = '00010';
            else
                res = '00011';
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

