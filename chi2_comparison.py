import sys, codecs, os, time
from scipy.stats import chi2_contingency
import numpy as np

os.environ['CLASSPATH'] = "../miralib/dist/*" # Need this to load the classes from the jar files

# in pyjnius 1.4+ we can use:
#import jnius_config
#jnius_config.add_options('-Xrs', '-Xmx4096')
#jnius_config.set_classpath('.', '../miralib/dist/*')

from jnius import autoclass
Log = autoclass("miralib.utils.Log")
Preferences = autoclass("miralib.utils.Preferences")
Project = autoclass("miralib.utils.Project")
DataRanges = autoclass("miralib.data.DataRanges")
DataSet = autoclass("miralib.data.DataSet")
Histogram = autoclass("miralib.data.ContingencyTable")
Similarity = autoclass("miralib.shannon.Similarity")
PValue = autoclass("miralib.shannon.PValue")

def getContingency(slice, proj):
    table = slice.getContingencyTable(proj)
    if not table.empty():
        rows = []
        for r in range(0, table.rowCount):
            row = table.getRow(r)
            if 0 in row: return None
            rows.append(row)
        return np.array(rows)
    else: 
        return None

Log.init()

# inputFile = "./data_nhanes/config.mira";
inputFile = "./data_lassa/config.mira";
outputFile = "./network.csv";

len = len(sys.argv)
if len == 2:
    inputFile = sys.argv[1];
else:
    for i in range(1, len):    
        if sys.argv[i] == "-in" and i + 1 < len:
            inputFile = sys.argv[i + 1];
        elif sys.argv[i] == "-out" and i + 1 < len:
            outputFile = sys.argv[i + 1];
        elif sys.argv[i] == "-miss" and i + 1 < len:
            preferences.missingString = sys.argv[i + 1];
        elif sys.argv[i] == "-mist" and i + 1 < len:
            preferences.missingThreshold = Project.stringToMissing(sys.argv[i + 1]); 
        elif sys.argv[i] == "-pval" and i + 1 < len:
            preferences.pValue = Project.stringToPValue(sys.argv[i + 1]);
        elif sys.argv[i] == "-algo" and i + 1 < len:
            preferences.depTest = Similarity.stringToAlgorithm(sys.argv[i + 1]);
        elif sys.argv[i] == "-surr" and i + 1 < len:
            preferences.surrCount = int(sys.argv[i + 1]);
        elif sys.argv[i] == "-chtr" and i + 1 < len:
            preferences.threshold = float(sys.argv[i + 1]);

preferences = Preferences()
project = Project(inputFile, preferences)
ranges = DataRanges();
data = DataSet(project);

print "Total number of data points:", data.getRowCount(ranges)
age = data.getVariableByAlias("Age at admission")
# age = data.getVariableByAlias("Age at Screening Adjudicated - Recode")
adults = age.createRange(18, 50)
ranges.update(age, adults)
print "Number of data points matching ranges:", data.getRowCount(ranges)

print "Total number of variables", data.getColumnCount()

data.removeColumns(data.getGroup("GWAS"))
data.removeColumns(data.getGroup("Treatment"))
data.removeColumns(data.getTable("Maximum Values of Lab Results"))
data.removeColumns(data.getTable("Last Day of Lab Results"))

# data.removeColumns(data.getGroup("Examination"))

count = data.getColumnCount()
output = [""] * (count + 1)
print "Number of selected variables", count

t0 = time.time()
print "*************************************"
total = 0
diff = 0
for i in range(0, count):
    vari = data.getColumn(i)
    if vari.categorical():
       print vari.getAlias()
       for j in range(0, count):
           varj = data.getColumn(j)
           if varj.categorical():
               slice = data.getSlice(vari, varj, ranges)
               if project.missingThreshold() <= slice.missing: continue
               obs = getContingency(slice, project)
               if obs == None: continue
               chi2, p, dof, ex = chi2_contingency(obs, correction=False)
               score = Similarity.calculate(slice, project.pvalue(), project)
               chi2_res = p < project.pvalue()
               simil_res = 0 < score
               total += 1
               if chi2_res != simil_res: 
                   mark = "<----"
                   diff += 1
               else: mark = ""
               print vari.getAlias(), varj.getAlias(), p, p < project.pvalue(), 0 < score, mark
t1 = time.time()

print "differences:", diff, "/" , total, "=",(100.0*diff/total), "%"
print "time:",t1-t0
