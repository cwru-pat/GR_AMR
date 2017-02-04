# Change these parameters
DB = "/home/chris/Desktop/blackhole-00040.visit/dumps.visit"
Mesh = "mesh"
Var = "DIFFK"
p0 = (46,50,50)
p1 = (54,50,50)
#########################################################################
OpenDatabase(DB)
AddPlot("Pseudocolor", Var)
TimeSliderSetState(20)
DrawPlots()
DefineScalarExpression("xc", "coord(%s)[0]" % Mesh)
Lineout(p0, p1, (Var))

# Get the data
SetActiveWindow(2)
SetActivePlots(0)
vals = GetPlotInformation()["Curve"]

f=open('/home/chris/Desktop/calib.txt','w')
# Write it as "x y z val"
for i in range(len(vals)/2+1 ):
  idx = i*2-1
  f.write("%g\n" % (vals[idx]))
  print "%g" % (vals[idx])
f.close()
