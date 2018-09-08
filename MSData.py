
class MSData:
    def __init__(self):
        self.fields = dict()
        self.prop = dict()
        self.fields['mz_list'] = []
        self.fields['intens_list'] = []
        self.fields['scan_number'] = []

    def size(self):
        return len(self.fields['intens_list'])

    def addMS(self,**attri):
        tmp = self.size()
        for feature in attri.keys():
            if not feature in self.fields.keys():
                self.fields[feature] = [float('nan') for i in range(tmp)]
            self.fields[feature].append(attri[feature])

    def setProp(self,**attri):
        for property in attri.keys():
            self.prop[property] = attri[property]

    def getProp(self,property):
        if not property in self.prop.keys():
            return None
        else:
            return self.prop[property]

    def getField(self,field):
        if field in self.fields.keys():
            return self.fields[field]
        else:
            return None

    def getSample(self,sn):
        I = self.fields['scan_number'].index(sn)
        if I < (self.size() - 1):
            tmp = []
            for feature in self.fields:
                tmp.append(self.fields[feature][I])
            return tmp
        else:
            return None

    def get(self,field,index):
        if field in self.fields.keys() and index < (self.size() - 1):
            return self.fields[field][index]
        else:
            return None



