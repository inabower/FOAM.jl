from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import numpy as np
import PyFoam

def readDict(fileName):
    #print('ParsedParameterFile(fileName).content')
    return parseForJl(ParsedParameterFile(fileName).content)

def readBoundaryDict(fileName):
    arr = ParsedParameterFile(fileName, boundaryDict=True).content
    ans = {}
    name = ""
    for f in arr:
        if name == "":
            name = f
        else:
            ans[name] = parseForJl(f)
            name = ""
    return ans

def parseForJl(masterDict):
    ans = {}
    for key in masterDict.keys():
        if type(masterDict[key]) is PyFoam.Basics.DataStructures.DictProxy:
            ans[key] = parseForJl(masterDict[key])
        elif type(masterDict[key]) is PyFoam.Basics.DataStructures.Dimension:
            ans[key] = masterDict[key].dims
        elif type(masterDict[key]) is PyFoam.Basics.DataStructures.Field:
            vals = masterDict[key]
            if vals.isUniform():
                valstr = str(masterDict[key]).replace("uniform","").strip()
                if valstr.find('(') > -1:
                    val = np.array([np.float(v) for v in valstr.replace("(","").replace(")","").split(' ')])
                else:
                    val = np.float(valstr)
                ans[key] = val
            else:
                ans[key] = np.array(masterDict[key])
        else:
            val = masterDict[key]
            try:
                val = float(val)
                if float(val).is_integer():
                    val = int(val)
            except:
                val = masterDict[key]
            ans[key] = val
    return ans
