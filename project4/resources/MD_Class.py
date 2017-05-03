import numpy as np
import random
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy.linalg as lg
from operator import mul
from Unit_Conversion import *


class MDSystem:
    structure = []
    timestep = 1.0
    dynamic = 0
    simulationTime = 0
    currentTime = 0
    dumpFilename = []
    dumpInterval = 1
    dumpMode = 0
    LogFilename = []
    LogInterval = 1
    LogMode = 0

    def __init__(self):
        # Dynamic: 1: NVE
        self.atoms = []

    def nve(self, timestep, simulationTime):

        self.timestep = cTime(timestep, 1)
        self.simulationTime = cTime(simulationTime, 1)
        while self.currentTime <= self.simulationTime:
            self.verlet()
            if (self.dumpMode != 0) & (np.mod(int(self.currentTime / self.timestep), self.dumpInterval) == 0):
                self.dump()
            if (self.LogMode != 0) & (np.mod(int(self.currentTime / self.timestep), self.LogInterval) == 0):
                self.Log()
            self.currentTime += self.timestep
    def force(self, distance):
        epsilon = 1  # 993.3   # 1.65 * 10^-21 J.                                    # eV
        sigma = 3.405  # A
        f = -24 * epsilon * (2 * (sigma ** 12 / distance ** 13) - (sigma ** 6 / distance ** 7))  # eV/A
        return f

    def PBC_distance(self, x1, x2):
        # Read cellLength from self.cellLength
        cellLength = self.structure.cellLength
        distance = x1 - x2
        # print 'distance = ', distance
        correctD = np.zeros((3))
        for i in [0, 1, 2]:
            correctD[i] = int(2 * (distance[i]) / cellLength[i]) * cellLength[i]
        distance =  correctD - distance
        # print 'correctD = ', correctD
        # print 'distance = ', distance
        return distance

    def checkPBC(self):
        cellLength = self.structure.cellLength
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            currentPosition = self.structure.position[i, :]
            for j in [0, 1, 2]:
                self.structure.position[i, j] -= math.floor(currentPosition[j] / cellLength[j]) * cellLength[j]

    def addStructure(self, newStructure):
        self.structure = newStructure

    def verlet(self):
        # Read positions and velocity
        # Apply verlet
        # Update position, velocity and time step
        # print 'verlet'
        self.updataForce()
        velocity2 = self.updataVelocity(self.structure.velocity, self.structure.force)
        self.updataPosition(velocity2)
        self.updataForce
        self.structure.velocity = self.updataVelocity(velocity2, self.structure.force)

    # def newton(self):
    #     h = self.timestep
    #     self.updataForce()
    #      for i in xrange(0, self.structure.totNumberOfAtom, 1):
    #          self.structure.position[i, :] +=


        # self.updataForce()
        # force1 = self.structure.force
        # self.updataPosition()
        # self.updataForce()
        # force2 = self.structure.force
        # self.updataVelocity(force1, force2)
        self.checkPBC()

    def setDump(self, filename, interval, mode):
        # mode = 1: position 2: position + velocity 3: position + velocity + force
        self.dumpFilename = filename
        self.dumpInterval = interval
        self.dumpMode = mode

    def setLog(self, filename, interval, mode):
        # mode = 1: position 2: position + velocity 3: position + velocity + force
        self.LogFilename = filename
        self.LogInterval = interval
        self.LogMode = mode
        f = open(self.LogFilename, 'a')
        f.write('time     ke      pe    etot  momentum diffCoeff temp \n')

    def dump(self):
        # write output file
        if self.dumpMode == 1:
            f = open(self.dumpFilename, 'a')
            f.write(str(self.structure.totNumberOfAtom))
            f.write('\n')
            f.write('Atom. timestep:')
            f.write(str(cTime(self.currentTime, 2)))
            f.write('\n')
            for i in xrange(0, self.structure.totNumberOfAtom, 1):
                f.write(str(self.structure.element[i][0]))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][0], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][1], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][2], 2)))
                f.write('\n')
        if self.dumpMode == 2:
            f = open(self.dumpFilename, 'a')
            f.write(str(self.structure.totNumberOfAtom))
            f.write('\n')
            f.write('Atom. timestep:')
            f.write(str(cTime(self.currentTime, 2)))
            f.write('\n')
            for i in xrange(0, self.structure.totNumberOfAtom, 1):
                f.write(str(self.structure.element[i][0]))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][0], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][1], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][2], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][0], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][1], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][2], 2)))
                f.write('\n')
        if self.dumpMode == 3:
            f = open(self.dumpFilename, 'a')
            f.write(str(self.structure.totNumberOfAtom))
            f.write('\n')
            f.write('Atom. timestep:')
            f.write(str(cTime(self.currentTime, 2)))
            f.write('\n')
            for i in xrange(0, self.structure.totNumberOfAtom, 1):
                f.write(str(self.structure.element[i][0]))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][0], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][1], 2)))
                f.write(' ')
                f.write(str(cLength(self.structure.position[i, :][2], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][0], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][1], 2)))
                f.write(' ')
                f.write(str(cVelocity(self.structure.velocity[i, :][2], 2)))
                f.write(' ')
                f.write(str(cForce(self.structure.force[i, :][0], 2)))
                f.write(' ')
                f.write(str(cForce(self.structure.force[i, :][1], 2)))
                f.write(' ')
                f.write(str(cForce(self.structure.force[i, :][2], 2)))
                f.write('\n')

        print 'dump', cTime(self.currentTime, 2)

    def Log(self):
        # Write log file
        Ek = cEnergy(self.kineticEnergy(), 2)
        Ep = cEnergy(self.potentialEnergy(), 2)
        f = open(self.LogFilename, 'a')
        f.write(str(cTime(self.currentTime, 2)))
        f.write(' ')
        f.write(str(Ek))
        f.write(' ')
        f.write(str(Ep))
        f.write(' ')
        f.write(str(Ek + Ep))
        f.write(' ')
        f.write(str(self.momentum()))
        f.write(' ')
        f.write(str(cDiffCoeff(self.DiffCoef(), 2)))
        f.write(' ')
        f.write(str(cTemperature(self.temperature(), 2)))
        f.write('\n')
        print 'Time = ', cTime(self.currentTime, 2), 'Ek = ', Ek, 'Ep = ', Ep, 'Etot = ', Ek + Ep, 'p = ', self.momentum(), 'DiffCoeff = ', cDiffCoeff(self.DiffCoef(), 2), 'T = ', cTemperature(self.temperature(), 2)

    def updataForce(self):
        # evaluate force from current position
        self.structure.force = np.zeros([self.structure.totNumberOfAtom, 3])
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            for j in xrange(i + 1, self.structure.totNumberOfAtom, 1):
                atomDistance = self.PBC_distance(self.structure.position[i, :], self.structure.position[j, :])
                atomDirection = atomDistance / lg.norm(atomDistance)
                atomForce = self.force(lg.norm(atomDistance))
                self.structure.force[i, :] += atomForce * atomDirection
                self.structure.force[j, :] += -atomForce * atomDirection
                # if i == 1:
                #     print 'i = ', i,'Position = ',self.structure.position[i]
                #     print 'j = ', j, 'Position = ',self.structure.position[j],
                #     print 'distance = ',atomDistance, 'normD = ', lg.norm(atomDistance)
                #     print 'atomForce = ', atomForce
                #     print 'force i = ', self.structure.force[i, :], 'force j = ', self.structure.force[j, :]



    def updataPosition(self, v):
        # updata position from current velocity and force
        h = self.timestep
        position_1 = self.structure.position
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            position_1[i, :] = self.structure.position[i, :] + h * v[i, :]

    def updataVelocity(self, v, f):
        v2 = np.zeros([self.structure.totNumberOfAtom, 3])
        h = self.timestep
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            v2[i, :] = v[i, :] + f[i, :] * h / 2 / self.structure.mass[i]
        return v2
        # def updataVelocity(self, force1, force2):
        #     velocity1 = self.structure.velocity
        #     h = self.timestep
        #     for i in xrange(0, self.structure.totNumberOfAtom, 1):
        #         self.structure.velocity[i, :] = velocity1[i, :] + h/2.0*(force1[i, :]+force2[i, :])/self.structure.mass[i]

    def kineticEnergy(self):
        totEk = 0
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            totEk += 0.5 * self.structure.mass[i] * lg.norm(self.structure.velocity[i]) ** 2
        return totEk

    def temperature(self):
        return 2.0 / 3.0 * self.kineticEnergy() / self.structure.totNumberOfAtom

    def potentialEnergy(self):
        totEp = 0
        epsilon = 1  # 993.3   # 1.65 * 10^-21 J.                                    # eV
        sigma = 3.405  # A
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            for j in xrange(i + 1, self.structure.totNumberOfAtom, 1):
                distance = lg.norm(self.PBC_distance(self.structure.position[i, :], self.structure.position[j, :]))
                totEp += 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)
        totEp /= 2
        return totEp

    def momentum(self):
        totMem = 0
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            totMem += sum(self.structure.velocity[i, :])
        return totMem

    def MSD(self):
        rt = 0.0
        for i in xrange(0, self.structure.totNumberOfAtom, 1):
            # print self.structure.position[i, :], self.structure.initStructure[i, :]
            rt += (lg.norm(self.PBC_distance(self.structure.position[i, :], self.structure.initStructure[i, :])))**2

        # print 'rt = ', rt
        rt /= self.structure.totNumberOfAtom
        return rt

    def DiffCoef(self):
        if(self.currentTime == 0):
            return 0
        else:
            return self.MSD()/6/self.currentTime


class Atom:
    cellLength = []  # A
    totNumberOfAtom = 0
    latConst = 0  # A
    position = np.array([0.0, 0.0, 0.0])  # A
    velocity = np.array([0.0, 0.0, 0.0])  # A/fs
    force = np.array([0.0, 0.0, 0.0])  # eV/A
    idx = np.array([0])
    mass = np.array([0.0])  # grams/mole
    element = np.array([0])
    initStructure = np.array([0.0, 0.0, 0.0])

    def addAtom(self, newIdx, newPosition, newVelocity, newMass, newElement):
        #  print self.idx
        #  print newIdx
        self.idx = np.vstack([self.idx, newIdx])
        self.position = np.vstack([self.position, newPosition])
        self.mass = np.vstack([self.mass, newMass])
        self.velocity = np.vstack([self.velocity, newVelocity])
        self.force = np.vstack([self.force, [0, 0, 0]])
        self.element = np.vstack([self.element, newElement])
        # print self.position
        # self.idx.append(newIdx)
        # self.position.append(newPosition)
        # self.mass.append(newMass)
        # self.velocity.append(newVelocity)

    def createFCC(self, element, mass):
        # CellLength adjustment
        self.supercell = np.divide(self.cellLength, self.latConst)
        self.supercell = self.supercell.astype(int)
        self.cellLength = self.supercell * self.latConst
        # Creating supercell
        FCCPosition = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]]) * self.latConst
        idx = 0
        self.totNumberOfAtom = reduce(mul, self.supercell) * 4
        for i in xrange(0, self.supercell[0], 1):
            for j in xrange(0, self.supercell[1], 1):
                for k in xrange(0, self.supercell[2], 1):
                    FCCPosition[:, 0] += i * self.latConst
                    FCCPosition[:, 1] += j * self.latConst
                    FCCPosition[:, 2] += k * self.latConst
                    Atom.addAtom(self, idx + 1, FCCPosition[0], [0, 0, 0], mass, element)
                    Atom.addAtom(self, idx + 2, FCCPosition[1], [0, 0, 0], mass, element)
                    Atom.addAtom(self, idx + 3, FCCPosition[2], [0, 0, 0], mass, element)
                    Atom.addAtom(self, idx + 4, FCCPosition[3], [0, 0, 0], mass, element)
                    idx += 4
                    FCCPosition = np.array([[0, 0, 0], [0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5]]) * self.latConst

    def createRand(self, numberOfAtom, element, mass):

        for i in xrange(0, numberOfAtom, 1):
            # print i
            # print [random.uniform(0, self.cellLength[1]), random.uniform(0, self.cellLength[2]), random.uniform(0, self.cellLength[3])]
            # print mass
            self.totNumberOfAtom += 1
            Atom.addAtom(self, i + 1, [random.uniform(0, self.cellLength[0]), random.uniform(0, self.cellLength[1]),
                                       random.uniform(0, self.cellLength[2])], [0.0, 0.0, 0.0], mass, element)

    def __init__(self, numberOfAtom, cellLength, mass, element, latConst, mode):
        # unit: cell length: A, mass: g/mol
        self.cellLength = cellLength
        if mode == 2:
            self.latConst = latConst
            Atom.createFCC(self, element, cMass(mass, 1))
        if mode == 1:
            Atom.createRand(self, numberOfAtom, element, cMass(mass, 1))
        self.idx = np.delete(self.idx, 0)
        self.mass = np.delete(self.mass, 0)
        self.position = np.delete(self.position, 0, axis=0)
        self.velocity = np.delete(self.velocity, 0, axis=0)
        self.force = np.delete(self.force, 0, axis=0)
        self.element = np.delete(self.element, 0, axis=0)
        self.initStructure = np.copy(self.position)
        # print N

    def printStructure(self):
        print self.velocity[1]
        print cVelocity(self.velocity[1], 2)
        print 'id mass x y z vx vy vz fx fy fz'
        for i in xrange(0, self.totNumberOfAtom, 1):
            print self.idx[i], cMass(self.mass[i], 2), cLength(self.position[i], 2), cVelocity(self.velocity[i],
                                                                                               2), cForce(self.force[i],
                                                                                                          2)

    def showStructures(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in xrange(0, self.totNumberOfAtom, 1):
            ax.scatter(cLength(self.position[i, 0], 2), cLength(self.position[i, 1], 2),
                       cLength(self.position[i, 2], 2))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()

    def assignTemperture(self, t):
        # print 'Assign temperature'
        t = cTemperature(t, 1)
        kEcurrent = 0
        kb = 1  # eV/K
        sigma = (kb * t / self.mass[0]) ** 0.5
        mu = 0
        totMem = 0
        for i in xrange(0, self.totNumberOfAtom, 1):
            nv = 0
            self.velocity[i, :] = np.random.normal(mu, sigma, 3)
            totMem += sum(self.velocity[i, :])
        print totMem
        # Remove momentum
        difference = totMem / (3 * self.totNumberOfAtom)
        for i in xrange(0, self.totNumberOfAtom, 1):
            self.velocity[i, :] = self.velocity[i, :] - [difference, difference, difference]

