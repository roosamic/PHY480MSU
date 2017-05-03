import numpy as np
from MD_Class import *

# myCell = Atom(10, [20, 20, 20], 39.948, 1, 5.260, 2)
myCell8 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell9 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell10 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell11 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell12 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell13 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell14 = Atom(1, [20, 20, 20], 39.948, 1, 5.405, 2)
myCell15 = Atom(1, [20, 20, 20], 39.948, 1, 5.260, 2)

myCell8.assignTemperture(300)
myCell9.assignTemperture(100)
myCell10.assignTemperture(50)
myCell11.assignTemperture(0.01)
# myCell12.assignTemperture(110)
# myCell13.assignTemperture(130)
# myCell14.assignTemperture(140)
# myCell15.assignTemperture(10)


s8 = MDSystem()
s8.addStructure(myCell8)
s8.setDump('traj_300.xyz', 100, 3)
s8.setLog('Ar_log_300.txt', 100, 1)

# s9 = MDSystem()
# s9.addStructure(myCell9)
# s9.setDump('traj9_100.xyz', 100, 3)
# s9.setLog('Ar_log9_100.txt', 100, 1)
#
# s10 = MDSystem()
# s10.addStructure(myCell10)
# s10.setDump('traj10_50.xyz', 100, 3)
# s10.setLog('Ar_log10_50.txt', 100, 1)
#
# s11 = MDSystem()
# s11.addStructure(myCell11)
# s11.setDump('traj11_0.xyz', 100, 3)
# s11.setLog('Ar_log11_0.txt', 100, 1)
#
# s12 = MDSystem()
# s12.addStructure(myCell12)
# s12.setDump('traj12_110.xyz', 100, 3)
# s12.setLog('Ar_log12_110.txt', 100, 1)
#
# s13 = MDSystem()
# s13.addStructure(myCell13)
# s13.setDump('traj13_130.xyz', 100, 3)
# s13.setLog('Ar_log13_130.txt', 100, 1)
#
# s14 = MDSystem()
# s14.addStructure(myCell14)
# s14.setDump('traj14_140.xyz', 100, 3)
# s14.setLog('Ar_log14_140.txt', 100, 1)
#
# s15 = MDSystem()
# s15.addStructure(myCell15)
# s15.setDump('traj15_005.xyz', 20, 3)
# s15.setLog('Ar_log15_005.txt', 20, 1)


s8.nve(0.001, 10)
# s9.nve(0.001, 50)
# s10.nve(0.001, 50)
# s11.nve(0.001, 50)
# s12.nve(0.001, 50)
# s13.nve(0.001, 50)
# s14.nve(0.001, 50)
# s15.nve(0.05, 50)



