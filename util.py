import re
from pymol import cmd
from pymol import stored
from pymol import selector
import subprocess
import numpy as np
from . import tempdir
import os
from sys import platform

viewer = 'primus'
if platform == 'win32' or platform == 'win64':
    viewer = 'primusqt'


# run crysol in predictive mode for a given selection
def predcrysol(modelName, crycalc):
    """Predicts scattering from the currently active model"""
    with tempdir.TemporaryDirectory() as tmpdir:
        if os.path.isfile(modelName + ".pdb"):
            print("Warning: PDB file already exists! Will be rewritten.")
        cmd.save(modelName + ".pdb", modelName)
        pdbFullPath = os.path.abspath(modelName + ".pdb")
        print("{} is saved to {}".format(modelName, pdbFullPath))
        if crycalc == "yes":
            print("CRYSOL calculation using explicit hydrogens")
            systemCommand(["crysol", "-eh", pdbFullPath])
        else:
            systemCommand(["crysol", pdbFullPath])
        result = parseCrysolLog(modelName + "00.log")
        Rg = result['Rg']
        eDens = result['eDens']
        df = tmpdir.move_out_numbered(modelName + "00.int", modelName, '.int')

    print("CRYSOL Theoretical Rg = " + repr(Rg))
    print("CRYSOL Average electron density = " + repr(eDens))
    print(".int file written to " + df)
    return df

    # run crysol in fit mode


def fitcrysol(modelName, dataName, crycalc, showFit):
    if dataName is None:
        print("Please set the SAXS .dat file to fit")
        return False
    if not os.path.isfile(dataName):
        print("SAXS .dat file \'" + dataName + "\' not found")
        return False
    if not os.path.isfile(modelName + ".pdb"):
        print("PDB file \'" + modelName + ".pdb" + "\' not found")
        return False
    fileFullPath = os.path.abspath(dataName)
    pdbFullPath = os.path.abspath(modelName + ".pdb")

    with tempdir.TemporaryDirectory() as tmpdir:
        if crycalc == "yes":
            print("CRYSOL calculation using explicit hydrogens")
            systemCommand(["crysol", "-eh", pdbFullPath, fileFullPath])
        else:
            systemCommand(["crysol", pdbFullPath, fileFullPath])
        logfile = modelName + "00.log"
        result = parseCrysolLog(logfile)
        Rg = result['Rg']
        chi2 = result['chi2']
        eDens = result['eDens']
        df = tmpdir.move_out_numbered(modelName + "00.fit", modelName, '.fit')
    fitResult = result
    fit = os.path.basename(df)
    # logfn = tmpdir.move_out_numbered(logfile, fid, '.log')

    print(".log file written to " + logfile)
    print(".fit file written to " + df)

    print("CRYSOL Theoretical Rg = " + repr(Rg))
    print("CRYSOL Chi-square = " + repr(chi2))
    print("CRYSOL Average electron density = " + repr(eDens))
    if showFit:
        systemCommand([viewer, fit])
    return fit, fitResult


def crysolRefinementSalipro(rot_min_ang, rot_max_ang, rot_step_ang,
                            scaffold_min, scaffold_max, scaffold_step,
                            protName, membName, scafName, dataName,
                            prefixName, runNumber):
    """Refine the membrane protein lipids scaffolding proteins
     complex against experimental data"""
    angs = np.arange(rot_min_ang, rot_max_ang, rot_step_ang)
    numScaffoldCopies = np.arange(scaffold_min, scaffold_max, scaffold_step)
    res = {"angle": -9999, "number-of-scaffolds": -9999, "chi2": 9999}
    with tempdir.TemporaryDirectory() as tmpdir:
        # copy data file
        tmpdir.copy_in(dataName)
        best = ""
        fitBest = ""
        for counter1, ang in enumerate(angs):
            for counter2, num in enumerate(numScaffoldCopies):
                cmd.refresh()
                modelName = builderSalipro(protName, scafName, membName, prefixName, runNumber, num, ang, True)
                runNumber += 1
                if modelName == "bad model":
                    print("Bad model parameters: ang: {} num: {}".format(ang, num))
                    continue
                cmd.save(modelName + ".pdb", modelName)
                fit, fitResult = fitcrysol(modelName, os.path.basename(dataName), "yes", False)
                cmd.wizard("message",
                           "Refinement: {} ".format(1 + counter2 + counter1 * (len(numScaffoldCopies))) +
                           " out of {} steps. Chi2: {}".format(len(numScaffoldCopies) * len(angs), fitResult['chi2']))
                if float(fitResult['chi2']) < float(res['chi2']):
                    if best != "": cmd.delete(best)
                    res['chi2'] = fitResult['chi2']
                    res['number-of-scaffolds'] = num
                    res['angle'] = ang
                    best = modelName
                    fitBest = fit
                else:
                    cmd.delete(modelName)
        tmpdir.move_out(best + ".pdb")
        tmpdir.move_out(fitBest)
    cmd.wizard()
    print("Best model: Number of Scaffolds = {}; Angle = {}".format(res['number-of-scaffolds'], res['angle']))
    print("Chi^2 : {} Best model name : {}".format(res['chi2'], best))
    return best, fitBest, counter1*counter2

    # run crysol in fit mode for detergents


def crysolRefinementDetergent(rot_min_ang, rot_max_ang, rot_step_ang,
                              dens_min_ang, dens_max_ang, dens_step_ang,
                              protName, membName, dataName, prefixName, runNumber):
    """Refine the membrane protein detergent complex against experimental data"""
    angs = np.arange(rot_min_ang, rot_max_ang, rot_step_ang)
    dens = np.arange(dens_min_ang, dens_max_ang, dens_step_ang)
    res = {"angle": -9999, "number-of-scaffolds": -9999, "chi2": 9999}
    with tempdir.TemporaryDirectory() as tmpdir:
        # copy data file
        tmpdir.copy_in(dataName)
        best = ""
        fitBest = ""

        for counter1, ang in enumerate(angs):
            for counter2, densAng in enumerate(dens):
                cmd.refresh()
                modelName = builderDetergent(protName, membName, prefixName, runNumber, ang, densAng, True)
                runNumber += 1
                if modelName == "bad model":
                    print("Bad model parameters: ang: {} densAng: {}".format(ang, densAng))
                    continue
                cmd.save(modelName + ".pdb", modelName)
                fit, fitResult = fitcrysol(modelName, os.path.basename(dataName), "yes", False)
                cmd.wizard("message",
                           "Refinement: {} ".format(1 + counter2 + counter1 * (len(dens))) +
                           " out of {} steps. Chi2: {}".format(len(dens) * len(angs), fitResult['chi2']))
                # if model fits better - store it
                if float(fitResult['chi2']) < float(res['chi2']):
                    if best != "": cmd.delete(best)
                    res['chi2'] = fitResult['chi2']
                    res['lipid density'] = densAng
                    res['max-polar-angle'] = ang
                    best = modelName
                    fitBest = fit
                else:
                    cmd.delete(modelName)
        if res['chi2'] < 9999:
            tmpdir.move_out(best + ".pdb")
            tmpdir.move_out(fitBest)
        else:
            return "Bad parameters", "No good fit found!"
    print("Best model: Lipid Density : {} Max Polar Angle = {})".format(res['lipid density'], res['max-polar-angle']))
    print("Chi^2 : {} Best model name : {}".format(res['chi2'], best))
    return best, fitBest, counter1*counter2

    # run crysol in fit mode for nanodisc


def crysolRefinementNanodisc(x_min, x_max, x_step, y_min, y_max, y_step,
                             protName, membName, scafName, dataName, prefixName, runNumber):
    """Refine the membrane protein detergent complex against experimental data"""
    xs = np.arange(x_min, x_max, x_step)
    ys = np.arange(y_min, y_max, y_step)
    res = {"x-offset": -9999, "y-offset": -9999, "chi2": 9999}
    with tempdir.TemporaryDirectory() as tmpdir:
        # copy data file
        tmpdir.copy_in(dataName)
        best = ""
        fitBest = ""

        for counter1, x in enumerate(xs):
            for counter2, y in enumerate(ys):
                cmd.refresh()
                modelName = builderNanodisc(protName, membName, scafName, prefixName, runNumber, x, y, True)
                runNumber += 1
                if modelName == "bad model":
                    # print("Bad model parameters: " + str(z))
                    continue
                cmd.save(modelName + ".pdb", modelName)
                fit, fitResult = fitcrysol(modelName, os.path.basename(dataName), "yes", False)
                cmd.wizard("message",
                           "Refinement: {} ".format(1 + counter2 + counter1 * (len(ys))) +
                           " out of {} steps. Chi2: {}".format(len(xs) * len(ys), fitResult['chi2']))
                # if model fits better - store it
                if float(fitResult['chi2']) < float(res['chi2']):
                    if best != "": cmd.delete(best)
                    res['chi2'] = fitResult['chi2']
                    res['x-offset'] = x
                    res['y-offset'] = y
                    best = modelName
                    fitBest = fit
                else:
                    cmd.delete(modelName)
        tmpdir.move_out(best + ".pdb")
        tmpdir.move_out(fitBest)
    print("Best model: Offset coordinates in XY plane : ({},{})".format(res['x-offset'], res['y-offset']))
    print("Chi^2 : {} Best model name : {}".format(res['chi2'], best))
    return best, fitBest, counter1*counter2


def builderSalipro(protein, scaffold, membrane, prefixName, runNumber, n_sym=9, initRotAngle=45, refine=False):
    """
    builds and refines MP-salipro systems
    """
    # Checking time of builder function execution
    if protein != None:
        print('protein is: ' + protein)
    print('scaffold is: ' + scaffold)
    print('membrane is: ' + membrane)
    print('copies of scaffold: ' + str(n_sym))

    print('rotation angle of scaffold is: ' + str(initRotAngle))

    empty = False
    if protein is None: empty = True
    print("State of empty/not-empty: {}".format(empty))
    cmd.reset()
    tmp_prot = "tmp_prot" + str(runNumber)
    tmp_scaffold = "tmp_scaffold" + str(runNumber)
    tmp_memb = "tmp_memb" + str(runNumber)
    tmp_origin = "origin" + str(runNumber)
    # copies to delete later
    cmd.copy(tmp_scaffold, scaffold)  # store initial
    cmd.copy(tmp_memb, membrane)  # store initial
    if not empty:
        cmd.copy(tmp_prot, protein)  # store initial
    center(tmp_memb)
    center(tmp_scaffold)
    cmd.pseudoatom(tmp_origin, pos=[0, 0, 0])
    cmd.origin(tmp_origin)
    t_sap = 10  # findMaxDist("tmp_scaffold")/2.0     # approximate half thickness of saposin monomer

    if not empty:
        avXY = TMdistCheck(tmp_prot, 0.2)
        if avXY == -1: return "bad model"
        radXY = avXY / 2.0
        # remove lipids inside pore
        cmd.remove("br. {} within {} of {}".format(tmp_memb, radXY, tmp_origin))
        r_lipHead = 4.7  # Area(POPC) = 65 A^2 => radius = sqrt(A/pi); reasonable estimate
        numLipLayers = 2.0  # number of lipid layers between TM of core and Saposin tmp_scaffold
        inRadius = avXY + numLipLayers * r_lipHead + t_sap  # equatorial position of tmp_scaffold center of mass
    else:
        inRadius = 27.0  # equatorial position of tmp_scaffold center of mass
    outRadius = inRadius + t_sap  # cut-lipids beyond this distance
    print("Inner radius: {}".format(inRadius))
    print("Outer radius: {}".format(outRadius))

    n_sym = int(n_sym)
    rotAng = 360. / float(n_sym)

    # chain ID
    id = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    s0 = range(0, n_sym, 1)
    # Build symmates with desired rotations
    cmd.rotate("x", "90", tmp_scaffold)
    for i in s0:
        angle = + i * rotAng
        cmd.copy("seg{}".format(i), tmp_scaffold)
        cmd.origin(tmp_origin)
        cmd.rotate("y", str(initRotAngle), "seg{}".format(i))
        cmd.origin(tmp_origin)
        cmd.translate("[0,{},0]".format(inRadius), "seg{}".format(i))
        chn = id[i + len(s0)]
        cmd.alter("seg{}".format(i), "chain = '{}'".format(chn))
        cmd.origin(tmp_origin)
        cmd.rotate("z", str(angle), "seg{}".format(i))

    # remove lipids beyond border encase by saposins
    cmd.remove("br. org and {} beyond {} of {}".format(tmp_memb, outRadius, tmp_origin))
    # remove lipids clashing with tmp_protein core and saposins
    if not empty:
        cmd.remove("br. org and {} within 0.35 of pol. and not hydro".format(tmp_memb))
    else:
        cmd.remove("br. org and {} within 0.3 of seg* and not hydro".format(tmp_memb))

    # Combine into a single PyMol object
    if empty:
        s = "empty_{}_{}_{}".format(membrane, scaffold, int(initRotAngle))
    else:
        s = "{}_{}_{}_{}".format(protein, membrane, scaffold, int(initRotAngle))
    if refine:
        s = "{}_{}_{}_{}_{}".format(protein, membrane, scaffold, int(initRotAngle), int(n_sym))
    if prefixName:
        s = "{}{}".format(prefixName, s)
    # cmd.create(s, protein, tmp_memb, "seg*")
    cmd.create(s, "({},{}, seg*)".format(protein, tmp_memb))
    cmd.save(s + ".pdb", s)

    cmd.delete(tmp_memb)
    cmd.delete(tmp_scaffold)
    cmd.delete(tmp_prot)
    cmd.delete("seg*")
    cmd.delete(tmp_origin)

    return s


def builderNanodisc(protein, membrane, scaffold, prefixName, runNumber, x=0, y=0, refine=False):
    """
    builds a MP-nanodisc systems
    scaffold in this case is a double belt of MSP
    """
    # Checking time of builder function execution
    if protein != None:
        print('protein is: ' + protein)
    print('scaffold is: ' + scaffold)
    print('membrane is: ' + membrane)
    empty = False
    if protein is None: empty = True
    if not empty:
        tmp_prot = "tmp_prot" + str(runNumber)
        cmd.copy(tmp_prot, protein)  # store initial
        cmd.translate("[{},{},0]".format(x, y), tmp_prot)
    print("State of empty/not-empty: {}".format(empty))
    # copies to delete later
    tmp_scaffold = "tmp_scaffold" + str(runNumber)
    tmp_memb = "tmp_memb" + str(runNumber)
    tmp_origin = "origin" + str(runNumber)
    cmd.copy(tmp_scaffold, scaffold)  # store initial
    cmd.copy(tmp_memb, membrane)  # store initial

    center(tmp_memb)
    center(tmp_scaffold)
    cmd.pseudoatom(tmp_origin, pos=[0, 0, 0])
    cmd.origin(tmp_origin)
    outRadius = findAverDist(tmp_scaffold)
    print("Max distance from origin to scaffold in xy plane: {}".format(outRadius))
    # remove lipids beyond border encased by MSP
    cmd.remove("br. org and {} beyond {} of {}".format(tmp_memb, outRadius, tmp_origin))

    # remove lipids clashing with tmp_protein core
    if not empty:
        avXY = TMdistCheck(tmp_prot, 0.2)
        if avXY == -1: return "bad model"
        minXY = avXY / 2.0
        # remove lipids inside pore
        cmd.remove("br. org and {} within {} of {}".format(tmp_memb, minXY, tmp_origin))
        print("Mean distance if TM cross-section in xy plane: {}".format(avXY))

    if empty:
        cmd.remove("br. org and {} within 0.4 of {} and not hydro".format(tmp_memb, tmp_scaffold))
        s = "empty_{}_{}".format(membrane, scaffold)
    else:
        cmd.remove("br. org and {} within 0.3 of pol. and not hydro".format(tmp_memb))
        s = "{}_{}_{}".format(protein, membrane, scaffold)
    if refine: s += "{}_{}".format(int(x), int(y))
    if prefixName:
        s = "{}{}".format(prefixName, s)
    cmd.create(s, "({},{}, {})".format(protein, tmp_scaffold, tmp_memb))
    cmd.save(s + ".pdb", s)

    cmd.delete(tmp_memb)
    cmd.delete(tmp_scaffold)
    cmd.delete(tmp_prot)
    cmd.delete(tmp_origin)
    return s


def builderDetergent(protein, detergent, prefixName, runNumber, ang=None, densAng=None, refine=False):
    """
    builds MP - detergent complex using a single detergent molecule
    """

    # Checking object names
    print('detergent is: ' + detergent)
    tmp_deter = "tmp_deter"+str(runNumber)
    tmp_prot  = "tmp_prot" + str(runNumber)
    tmp_origin = "origin0" + str(runNumber)
    cmd.copy(tmp_deter, detergent)  # store initial
    # molecules initially aligned along Z-axis on import, need to rotate into XY plane for protocol: 
    print("Rotating initial {} aligned along Z ==> into XY plane...".format(detergent))
    cmd.rotate("x", -90, tmp_deter)
    center(tmp_deter)
    # if it is an empty assembly --> build a micelle
    if protein is None:
        # table with common empty micelle parameters
        # https://www.molecularbiophysics.physik.lmu.de/publications/Lipfert_etal_JPCB07.pdf
        if refine:
            radius = ang
            numberOfDetergents = densAng
        else:
            radius = findMaxDist(tmp_deter)
            numberOfDetergents = 300
        s = builderMicelle(tmp_deter, 2*radius, numberOfDetergents)
        cmd.save(s + ".pdb", s)
        return s

    print('protein   is: ' + protein)
    cmd.copy(tmp_prot, protein)  # store initial
    # Determine max distance of TM cross-section (xy plane)
    # r        = TMdistCheck("tmp_prot", 2.0)
    # if r == -1: return "bad model"
    detR = findMaxDist(tmp_deter)
    # print("Max distance if TM cross-section is in a xy plane: " + (str)(r))
    print("Max distance of detergent : " + str(detR))

    # Shrink detergent along z axis to match ry
    # stretch = r/detR
    # affineStretch(detergent, stretch)
    # find new
    # detR     = detR * stretch
    # print("Max distance of detergent after shrinking: " + (str)(detR))
    # Create a ring of detergents using spherical coordinates
    # TODO: find automatically?
    if refine:
        # stochastic
        # theta = random.sample(range(-ang, ang), 10)
        # phi = random.sample(range(0, 360), densAng)
        # geometrical
        theta = np.arange(-ang, ang, 3)
        phi = np.arange(0, 361, densAng)
    else:
        # stochastic
        # theta = random.sample(range(-20, 20), 10)
        # phi = random.sample(range(0, 361), 100)
        # geometrical
        theta = range(-14, 14, 3)
        phi = range(0, 361, 10)  # find angular step from average density?
    # builderCorona(theta, phi, "tmp_deter", r, detR)
    builderCorona(theta, phi, tmp_deter, tmp_prot, detR)

    # combine components into single PYMOL object
    if refine:
        s = "{}_{}_{:d}_{:d}".format(protein, detergent, int(ang), int(densAng))
    else:
        s = "{}_complex_with_{}".format(protein, detergent)
    #print(f"REfINE={refine} for {s}")
    if prefixName:
        s = "{}{}".format(prefixName, s)
    # affineStretch("corona", 1.1)
    cmd.create(s, "({}, corona)".format(protein))
    cmd.delete("corona")
    cmd.delete(tmp_prot)
    cmd.delete(tmp_deter)
    cmd.delete(tmp_origin)
    cmd.save(s + ".pdb", s)
    return s

def fibonacci_sphere(samples):
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points.append([x, y, z])

    return points


def builderMicelle(detergent, r, numberOfDetergents):
    refresh()
    i = 0
    numberOfDetergents = int(numberOfDetergents)
    # FIXME: if number of detergents > 360, molecules may clash in space
    points = fibonacci_sphere(numberOfDetergents)
    print("Rotating detergent molecule for build protocol...")
    for x,y,z in points:
        i += 1
        t = np.degrees(np.arccos(z))
        p = np.degrees(np.arctan(y/x))
        cmd.copy("seg{}".format(i), detergent)
        cmd.alter("seg{}".format(i), "resi={}".format(i))  # assign residue numbers
        # put to to knot of fibonacci grid
        center("seg{}".format(i))
        cmd.translate("[0,{},0]".format(r), "seg{}".format(i))
        if x < 0: cmd.rotate("z", 180, "seg{}".format(i))
        cmd.rotate("z", str(t), "seg{}".format(i))
        cmd.rotate("y", str(p), "seg{}".format(i))
        #cmd.translate("[{},{},{}]".format(r * x, r * y, r * z), "seg{}".format(i))
    # if numberOfDetergents > 360:
    #     x1 = range(-180, 180)
    #     x2 = range(0, 360)
    #     theta = ([random.choice(x1) for _ in range(numberOfDetergents)])
    #     phi = ([random.choice(x2) for _ in range(numberOfDetergents)])
    # else:
    #     theta = random.sample(range(-180, 180), numberOfDetergents)
    #     phi = random.sample(range(0, 360), numberOfDetergents)
    # for t, p in zip(theta, phi):
    #     i += 1
    #     cmd.copy("seg{}".format(i), detergent)
    #     cmd.alter("seg{}".format(i), "resi={}".format(i))  # assign residue numbers
    #     # randomly sample on a sphere
    #     cmd.translate("[0,{},0]".format(r), "seg{}".format(i))
    #     cmd.rotate("x", str(t), "seg{}".format(i))
    #     cmd.rotate("z", str(p), "seg{}".format(i))
    s = "micelle_{}_{}_{}".format(detergent, int(r), int(numberOfDetergents))
    cmd.create(s, "seg*")
    cmd.delete("seg*")
    # center(f"{s}")
    # cmd.show_as("sticks","org")
    # could be streched if necessary
    # affineStretch(s, 10)
    return s


def builderCorona(theta, fi, detergent, protein, detR):
    # Build symmates with desired rotations
    refresh()
    cmd.pseudoatom("origin0"+protein, pos=[0, 0, 0])
    thetaSteps = len(theta)
    angleVer = np.linspace(-90, 90, thetaSteps)
    i = 0
    roffi = []
    # find surface grid todo: choose better one of these two approaches
    # for f in fi:
    # find distances to surface
    # cmd.rotate("z", (str)(-f), protein)
    # xLine = "{} and z > {} and z < {}  and y > {} and y < {} and x > 0".format(protein, -detR, detR, -5, 5)
    # atoms = cmd.index(xLine)
    # dlist = []
    # if len(atoms) == 0:
    #     print("No atoms at azimuthal angle: Phi = {}.".format(f))
    #     continue
    #     # find R in only one cross-section
    # for at1 in cmd.index("origin0"):
    #     for at2 in atoms[::10]:
    #         dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
    #         dlist.append(dist)
    # r = max(dlist) # * np.cos(np.deg2rad(t))
    # roffi.append(r)
    # cmd.rotate("z", (str)(f), protein)
    r = 0
    cmd.rotate("z", -90, detergent)
    for t, a in zip(theta, angleVer):
        for n, f in enumerate(fi):
            ###################EXPERIMENTAL###################
            cmd.rotate("z", str(-f), protein)
            cmd.translate("[0,0,{}]".format(r * np.sin(np.deg2rad(t))), protein)
            xLine = "{} and z > {} and z < {}  and y > {} and y < {} and x > 0 and name CA".format(protein, -detR, detR,
                                                                                                   -1, 1)
            atoms = cmd.index(xLine)
            dlist = []
            if len(atoms) == 0:
                print("No atoms at azimuthal angle: Phi = {}.".format(f))
                continue
                # find R in only one cross-section
            for at1 in cmd.index("origin0"):
                for at2 in atoms:
                    dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
                    dlist.append(dist)
            r = max(dlist)  # * np.cos(np.deg2rad(t))
            roffi.append(r)
            cmd.rotate("z", str(f), protein)
            cmd.translate("[0,0,{}]".format(-r * np.sin(np.deg2rad(t))), protein)
            ####################################################
            # r = roffi[n] / np.cos(np.deg2rad(t))
            i += 1
            cmd.copy("seg{}".format(i), detergent)
            cmd.alter("seg{}".format(i), "resi={}".format(i))  # assign residue numbers
            # corona
            cmd.rotate("y", str(-a), "seg{}".format(i))
            cmd.translate("[{},0,0]".format(r + 0.6 * detR * np.cos(np.deg2rad(a))), "seg{}".format(i))
            cmd.translate("[0,0,{}]".format((r + detR) * np.sin(np.deg2rad(t))), "seg{}".format(i))
            cmd.rotate("z", str(f), "seg{}".format(i))
            # print(f"seg{i} phi = {f} theta = {t} Distance: {r}") #DEBUG

    cmd.create("corona", "seg*")
    cmd.delete("seg*")
    cmd.delete("origin0"+protein)


def builderMembrane(lipid, runNumber):
    """
    build membrane bilayer from single lipid PDB file
    """
    refresh()
    cmd.load(lipid + ".pdb", "start_lipid")
    cmd.alter("start_lipid", "chain = 'X'")
    cmd.alter("start_lipid", "segi = 'mema'")
#    cmd.rotate('x', 90, "start_lipid")
    dmax = findMaxDist("start_lipid")

    # create lipid copies and translate them to new position
    nlip = 20  # number of lipids forming edge of bilayer

    s0 = range(1, nlip, 1)
    s1 = range(1, nlip + 1, 1)  # excludes first lipid

    step_x = 0  # translation in x (TODO: automatic determination of spacing without clashes)
    step_y = 7
    step_z = 0
    step_x2 = 7
    step_y2 = 0
    step_z2 = 0

    for i in s1:
        # first column
        cmd.copy("lip{}".format(i), "start_lipid")  # row of lipids
        cmd.alter("lip{}".format(i), "resi={}".format(i))  # change residue numbers
        y = i * step_y
        cmd.translate("[{},{},{}]".format(step_x, y, step_z), "lip{}".format(i))
        # generate remaining rows/columns in same leaflet
        for j in s0:
            k = int(nlip) * i + j  # TODO: general counter to write correct lipid number
            cmd.copy("lip{}".format(k), "lip{}".format(i))  # adjacent row of lipids
            cmd.alter("lip{}".format(k), "resi={}".format(k))  # change residue numbers
            x2 = j * step_x2
            cmd.translate("[{},{},{}]".format(x2, step_y2, step_z2), "lip{}".format(k))
        cmd.sort()  # sort atom order
    # create second leaflet
    # simple method by creating a single leaflet object:
    cmd.create("mema", "(lip*)")
    cmd.delete("lip*")

    cmd.copy("memb", "mema")
    cmd.alter("memb", "segi = 'memb'")
    cmd.rotate("x", 180, "memb")
    cmd.translate("[0,0,{}]".format((-1.0 * (dmax + 0.5))), "memb")
    # cmd.color("yellow", "segi = 'mema'")
    # cmd.color("blue", "segi = 'memb'")
    cmd.translate("[3.5,3.5,0]", "memb")  # optional shift of single leaflet to avoid aliphatic clashes
    s = "{}_bilayer".format(lipid)
    cmd.create(s, "(mema,memb)")
    cmd.delete("mema ,memb, start_lipid")
    center(s)
    cmd.save(s + ".pdb", s)
    cmd.reset()
    return s


def builderBicelle(protein, membrane, detergent, prefixName, runNumber, refine=False, ang=None, densAng=None):
    """
    builds MP - bicelle complex using a membrane and a single detergent molecule
    """
    # Checking object names
    print('protein   is: ' + protein)
    print('membrane  is: ' + membrane)
    print('detergent is: ' + detergent)
    empty = False
    if protein is None: empty = True
    if not empty:
        tmp_prot = "tmp_prot" + str(runNumber)
        cmd.copy(tmp_prot, protein)  # store initial
    print("State of empty/not-empty: {}".format(empty))
    # copies to delete later
    tmp_memb = "tmp_memb" + str(runNumber)
    tmp_deter =  "tmp_deter"  + str(runNumber)
    cmd.copy(tmp_memb, membrane)  # store initial
    cmd.copy(tmp_deter, detergent)  # store initial

    center(tmp_memb)
    center(tmp_deter)
    tmp_origin = "origin" + str(runNumber)
    cmd.pseudoatom(tmp_origin, pos=[0, 0, 0])
    cmd.origin(tmp_origin)
    # Determine max distance of TM cross-section (xy plane)
    r = TMdistCheck(tmp_prot, 2.0)
    if r == -1: return "bad model"
    detR = findMaxDist(tmp_deter)
    print("Max distance if TM cross-section is in a xy plane: " + str(r))
    print("Max distance of detergent : " + str(detR))

    ## Shrink detergent along z axis to match ry
    # stretch = r / detR
    # affineStretch("tmp_deter", stretch)
    # find new
    # detR = detR * stretch
    # print(f"Max distance of detergent after shrinking: {detR}")
    # Create a ring of detergents using spherical coordinates
    # FIXME: find automatically?
    if refine:
        # stochastic
        # theta = random.sample(range(-ang, ang), 10)
        # phi = random.sample(range(0, 360), densAng)
        # geometrical
        theta = np.arange(-ang, ang, 3)
        phi = np.arange(0, 361, densAng)
    else:
        # stochastic
        # theta = random.sample(range(-20, 20), 10)
        # phi = random.sample(range(0, 361), 100)
        # geometrical
        theta = range(-14, 14, 3)
        phi = range(0, 361, 10)  # find angular step from average density?

    builderCorona(theta, phi, tmp_deter, r, detR)
    #affineStretch("corona", 1.1)
    # remove lipids inside pore
    cmd.remove("br. {} within {} of {}".format(tmp_memb, r / 2.0, tmp_origin))
    cmd.origin(tmp_origin)
    # remove lipids beyond border encased by MSP
    #print("org and tmp_memb beyond {} of origin0".format(r))
    cmd.remove("br. org and {} beyond {} of {}".format(tmp_memb, r, tmp_origin))
    # remove lipids clashing with tmp_protein core and MSP scaffold and combine into a single PyMol object
    cmd.remove("br. org and {} within 0.3 of pol. and not hydro".format(tmp_memb))
    s = "{}_{}_{}".format(protein, membrane, detergent)
    if prefixName: s = "{}{}".format(prefixName, s)
    cmd.create(s, "({}, corona, {})".format(protein, tmp_memb))
    cmd.save(s + ".pdb", s)
    cmd.delete("{}, {}, {}, corona".format(tmp_prot, tmp_memb, tmp_deter))
    cmd.delete(tmp_origin)
    return s


def TMdistCheck(selection, z):
    """
    determine cross section distances about origin of trans-membrane region of protein PDB
    (distance from origin to protein atoms in xy plane)
    """
    cmd.pseudoatom("origin0", pos=[0, 0, 10])
    dmax = []
    for ang in range(0, 181, 10):
        dlist = []
        cmd.rotate('z', ang, selection)
        xLine = "{} and z > {} and z < {}  and y > {} and y < {}".format(selection, -z, z, -z, z)
        atoms = cmd.index(xLine)
        if len(atoms) == 0:
            print("No atoms in {}.".format(xLine))
            continue
        # find R in only one cross-section
        for at1 in cmd.index("origin0"):
            for at2 in atoms:
                dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
                dlist.append(dist)
        dmax.append(max(dlist))
        cmd.rotate('z', -ang, selection)
    if len(dmax) == 0:
        meanXY = -1
    else:
        meanXY = np.mean(dmax)
    cmd.pseudoatom("origin0", pos=[0, 0, 0])
    return meanXY


def refresh():
    cmd.reset()
    cmd.show_as("cartoon", "pol")
    # cmd.show_as("sticks", "org")
    cmd.show_as("spheres", "org")


def findMaxDist(selection):
    """
    finds the longest distance within a molecule
    """
    dlist = []
    for at1 in cmd.index(selection)[::10]:
        for at2 in cmd.index(selection)[::10]:
            dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
            dlist.append(dist)
    dmax = max(dlist)
    return dmax


def findAverDist(selection):
    dlist = []
    cmd.pseudoatom("origin0", pos=[0, 0, 0])
    for at1 in cmd.index("origin0"):
        for at2 in cmd.index(selection)[::30]:
            dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
            dlist.append(dist)
    return np.mean(dlist)


# parse crysol log file
def parseCrysolLog(logFileName):
    """Parse Crysol log file, obtain Chi2, Rg and eDens"""
    # will not parse crysol_summary.txt, but the .log file
    # created for each individual run

    chi2 = 9999
    Rg = -9999
    eDens = -9999

    position = -1
    counter = 0
    with open(logFileName, 'r') as rf:
        for line in rf:
            counter += 1
            if re.match("(.*)Fitting parameters(.*)", line):
                print("line number: " + repr(counter))
                position = counter + 2
            if counter == position:
                if line[66:73] != "*******":
                    chi2 = float(line[66:73])
            if re.match("(.*)Rg from the slope of net intensity(.*)", line):
                Rg = float(line[59:65])
            if re.match("(.*)Average electron density(.*)", line):
                eDens = float(line[59:66])
    rf.close()
    return {'chi2': chi2, 'Rg': Rg, 'eDens': eDens}


def writePdb(sel, prefix=""):
    pdbfn = prefix + sel + ".pdb"
    npdbfn = pdbfn.replace(" or ", "")
    npdbfn = npdbfn.replace(" and ", "")

    try:
        npdbfn = npdbfn.translate(npdbfn.maketrans('', '', " "))
    except ImportError:
        pass
    else:
        try:
            npdbfn = npdbfn.translate(" ")
        except ImportError:
            pass
    print("{} is saved to {}".format(sel, npdbfn))
    cmd.save(npdbfn, sel)
    return npdbfn


def systemCommand(command, **kwargs):
    status = subprocess.call(command, **kwargs)
    if 0 != status:
        print("WARNING, something went wrong while executing:\n"
              + ' '.join(command))
    return status


def center(selection):
    # center a molecule
    com = cmd.centerofmass(selection, -1.0)
    shift = np.array2string(np.negative(np.array(com)), separator=', ')
    cmd.translate(shift, selection)
    print('Shifting {} to center'.format(selection))


cmd.extend("center", center)


def affineStretch(selection, stretch):
    # stretch molecule using affine transformations
    stored.altered = []
    cmd.iterate_state(1, selector.process(selection), "stored.altered.append([x,{}*y,z])".format(stretch))
    cmd.alter_state(1, selection, "(x,y,z) = stored.altered.pop(0)")


cmd.extend("affineStretch", affineStretch)


def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec


def loadModel(filename, pymolObject):
    if filename:
        try:
            cmd.delete(pymolObject)
            cmd.load(filename + ".pdb")
        except Exception as e:
            print("For file {} and {}: {}".format(filename, pymolObject, e))


# Not used at the moment!!
def affine(selection, sh_x=0.1, sh_y=0.1):
    """ 
    Perform distortion from circular geometry on scaffold
    Peform rotation-translation-rotation on a PDB
    1. load scaffold
    2. apply transfromation
    sh_x,sh_y = shear factors (x,y)
    """
    # set-up transformation

    shift0 = 0.0
    shift1 = 0.0
    shift2 = 0.0
    #    sh_x = 0.1 # shear factor for transformation in x
    #    sh_y = 0.1 # shear factor for transformation in y

    rotMat = [1.0, sh_x, 0.0, shift0,
              sh_y, 1.0, 0.0, shift1,
              0.0, 0.0, 1.0, shift2,
              0.0, 0.0, 0.0, 1.0]

    # perform the transformation
    return cmd.transform_selection(selection, rotMat)


def curveBilayer(radius=15.0, shift=1.0):
    """ 
    Raise/lower height of central bilayer lipids to create a curved bilayer surface
    """
    shiftTop = shift  # shift in + Z
    shiftBot = -1.0 * float(shift)  # shift in -Z
    cmd.translate([0.0, 0.0, shiftTop], "br. segi 'mema' within %s of origin0" % radius)  # Adjust height of bilayer
    cmd.translate([0.0, 0.0, shiftBot], "br. segi 'memb' within %s of origin0" % radius)  # Adjust height of bilayer
    print("Adjusted vertical position of central lipids by %s A" % shift)
