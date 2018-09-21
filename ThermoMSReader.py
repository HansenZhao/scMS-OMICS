import MSFileReader
import MSData
from comtypes import COMError

def read(file_name,method = 0):
    rawfile = MSFileReader.ThermoRawfile(file_name)
    scan_from = rawfile.GetFirstSpectrumNumber()
    scan_to = rawfile.GetLastSpectrumNumber()
    time_from = rawfile.GetStartTime()
    time_to = rawfile.GetEndTime()
    mass_from = rawfile.GetLowMass()
    mass_to = rawfile.GetHighMass()
    scan_number = scan_to - scan_from + 1
    inst_name = rawfile.GetInstName()
    isUser, tol_value, tol_unit = rawfile.GetMassTolerance()
    if tol_unit == 1:
        tol_unit = 'ppm'
    elif tol_unit == 0:
        tol_unit = 'mmu'
    elif tol_unit == 2:
        tol_unit = 'amu'

    tolerence = '{:.4f} {}'.format(tol_value,tol_unit)

    obj = MSData.MSData()
    obj.setProp(scan_range = [scan_from,scan_to],time_range=[time_from,time_to],mass_range = [mass_from,mass_to],
                scan_number = scan_number, inst_name = inst_name, tolerence = tolerence)

    if method == 0:
        obj.setProp(read_method = 'Peaks')
        for i in range(scan_from,scan_to+1):
            tmp = rawfile.GetMassListFromScanNum(i)
            if not tmp[0][0]:
                obj.addMS(mz_list=(), intens_list=(), ppm=(), resolution=(), filter_info=(),
                          scan_number=i, scan_time=rawfile.RTFromScanNum(i))
            else:
                try:
                    intens, mz, MMU, PPM, resolution = rawfile.GetMassPrecisionEstimate(i)
                except COMError:
                    print('failed to read: scan {}'.format(i))
                    obj.setProp(scan_number = obj.getProp('scan_number') - 1)
                    continue
                filter_info = rawfile.GetFilterForScanNum(i)
                obj.addMS(mz_list=mz,intens_list=intens,ppm = PPM,resolution=resolution,filter_info=filter_info,scan_number = i,
                          scan_time = rawfile.RTFromScanNum(i))
    elif method == 1:
        obj.setProp(read_method='Centronized')
        for i in range(scan_from, scan_to + 1):
            mass, tmp = rawfile.GetMassListFromScanNum(i,centroidResult=True)
            filter_info = rawfile.GetFilterForScanNum(i)
            obj.addMS(mz_list=mass[0], intens_list=mass[1], filter_info=filter_info,scan_number = i, scan_time=rawfile.RTFromScanNum(i))
    else:
        obj.setProp(read_method='Profile')
        for i in range(scan_from, scan_to + 1):
            mass, tmp = rawfile.GetMassListFromScanNum(i)
            filter_info = rawfile.GetFilterForScanNum(i)
            obj.addMS(mz_list=mass[0], intens_list=mass[1], filter_info=filter_info,scan_number = i, scan_time=rawfile.RTFromScanNum(i))

    return obj