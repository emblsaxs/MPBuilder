import re
from pymol import cmd
from pymol import stored
from pymol import selector
import subprocess
import numpy as np
from . import tempdir
import os
import random


# run crysol in predictive mode for a given selection
def predcrysol(modelName, crycalc):
    '''Predicts scattering from the currently active model'''
    with tempdir.TemporaryDirectory() as tmpdir:
        if True == os.path.isfile(modelName + ".pdb"):
            print(f"Warning: PDB file already exists! Will be rewritten.")
        cmd.save(modelName + ".pdb", modelName)
        pdbFullPath = os.path.abspath(modelName + ".pdb")
        print(f"{modelName} is saved to {pdbFullPath}")
        if (crycalc == "yes"):
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
def fitcrysol(modelName, dataName, crycalc,  showFit):
    if False == os.path.isfile(dataName):
        print("SAXS .dat file \'" + dataName + "\' not found")
        return False
    if False == os.path.isfile(modelName + ".pdb"):
        print("PDB file \'" + modelName + ".pdb" + "\' not found")
        return False
    fileFullPath = os.path.abspath(dataName)
    pdbFullPath = os.path.abspath(modelName + ".pdb")

    with tempdir.TemporaryDirectory() as tmpdir:
        if (crycalc == "yes"):
            print("CRYSOL calculation using explicit hydrogens")
            systemCommand(["crysol", "-eh", pdbFullPath, fileFullPath])
        else:
            systemCommand(["crysol",pdbFullPath, fileFullPath])
        logfile = modelName + "00.log"
        result = parseCrysolLog(logfile)
        Rg = result['Rg']
        chi2 = result['chi2']
        eDens = result['eDens']
        df = tmpdir.move_out_numbered(modelName + "00.fit", modelName, '.fit')
    fitResult = result
    fit       = os.path.basename(df)
    #logfn = tmpdir.move_out_numbered(logfile, fid, '.log')

    print(".log file written to " + logfile)
    print(".fit file written to " + df)

    print("CRYSOL Theoretical Rg = " + repr(Rg))
    print("CRYSOL Chi-square = " + repr(chi2))
    print("CRYSOL Average electron density = " + repr(eDens))
    if showFit:
        systemCommand(["primus", fit])
    return fit, fitResult

def crysolRefinementSalipro(rot_min_ang, rot_max_ang, rot_step_ang, \
                            scaffold_min, scaffold_max, scaffold_step, \
                            protName, membName, scafName,  dataName, \
                            prefixName):
    '''Refine the membrane protein lipids scaffolding proteins
     complex against experimental data'''
    angs = np.arange(rot_min_ang, rot_max_ang, rot_step_ang)
    numScaffoldCopies = np.arange(scaffold_min, scaffold_max, scaffold_step)
    res = {"angle": -9999, "number-of-scaffolds": -9999, "chi2": 9999}
    with tempdir.TemporaryDirectory() as tmpdir:
        # copy data file
        tmpdir.copy_in(dataName)
        best    = ""
        fitBest = ""
        for counter1, ang in enumerate(angs):
            for counter2, num in enumerate(numScaffoldCopies):
                cmd.refresh()
                modelName = builderSalipro(protName, scafName, membName, prefixName, num, ang, True)
                cmd.save(modelName + ".pdb", modelName)
                fit, fitResult = fitcrysol(modelName, dataName, "yes", False)
                cmd.wizard("message",
                           f"Refinement: {counter2 + counter1 * (len(numScaffoldCopies))} "
                           f"out of {len(numScaffoldCopies) * len(angs)} steps. Chi2: {fitResult['chi2']}")
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
    print(f"Best model: Number of Scaffolds = {res['number-of-scaffolds']}; Angle = {res['angle']}")
    print(f"Chi^2 : {res['chi2']} Best model name : {best}")
    return best, fitBest


    # run crysol in fit mode for detergents
def crysolRefinementDetergent(rot_min_ang, rot_max_ang, rot_step_ang, \
                              dens_min_ang, dens_max_ang, dens_step_ang, \
                              protName, membName, dataName, prefixName):
    '''Refine the membrane protein detergent complex against experimental data'''
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
                modelName = builderDetergent(protName, membName, prefixName, ang, densAng, True)
                cmd.save(modelName + ".pdb", modelName)
                fit, fitResult = fitcrysol(modelName, dataName, "yes", False)
                cmd.wizard("message",
                           f"Refinement: {counter2 + counter1 * (len(dens))} out of {len(dens) * len(angs)} steps."
                           f"Chi2: {fitResult['chi2']}")
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
        tmpdir.move_out(best + ".pdb")
        tmpdir.move_out(fitBest)
    print(f"Best model: Lipid Density : {res['lipid density']} Max Polar Angle = {res['max-polar-angle']})")
    print(f"Chi^2 : {res['chi2']} Best model name : {best}")
    return best, fitBest

    # run crysol in fit mode for nanodisc
def crysolRefinementNanodisc(z_min, z_max, z_step, \
                              protName, membName, scafName, dataName, prefixName):
    '''Refine the membrane protein detergent complex against experimental data'''
    zs = np.arange(z_min, z_max, z_step)
    res = {"angle": -9999, "number-of-scaffolds": -9999, "chi2": 9999}
    with tempdir.TemporaryDirectory() as tmpdir:
        # copy data file
        tmpdir.copy_in(dataName)
        best = ""
        fitBest = ""

        for counter1, z in enumerate(zs):
            cmd.refresh()
            modelName = builderNanodisc(protName, membName, scafName, prefixName, z, True)
            cmd.save(modelName + ".pdb", modelName)
            fit, fitResult = fitcrysol(modelName, dataName, "yes", False)
            cmd.wizard("message",
                       f"Refinement: {counter1 } out of {len(zs)} steps."
                       f"Chi2: {fitResult['chi2']}")
            # if model fits better - store it
            if float(fitResult['chi2']) < float(res['chi2']):
                if best != "": cmd.delete(best)
                res['chi2'] = fitResult['chi2']
                res['vertical-offset'] = z
                best = modelName
                fitBest = fit
            else:
                cmd.delete(modelName)
        tmpdir.move_out(best + ".pdb")
        tmpdir.move_out(fitBest)
    print(f"Best model: Vertical Offset : {res['vertical-offset']}")
    print(f"Chi^2 : {res['chi2']} Best model name : {best}")
    return best, fitBest

def builderSalipro(protein, scaffold, membrane, prefixName, n_sym=9, initRotAngle=45, refine = False):
    """
    builds and refines MP-salipro systems
    """
    # Checking time of builder function execution
    print(f'protein is: {protein}')
    print(f'scaffold is: {scaffold}')
    print(f'membrane is: {membrane}')
    print(f'copies of scaffold: {n_sym}')

    print(f'rotation angle of scaffold is: {initRotAngle}')

    empty = False
    if protein == None: empty = True
    cmd.reset()
    # copies to delete later
    cmd.copy("tmp_scaffold",scaffold) # store initial
    cmd.copy("tmp_memb",membrane) # store initial
    cmd.copy("tmp_prot",protein) # store initial
    center("tmp_memb")
    center("tmp_scaffold")
    cmd.pseudoatom("origin0", pos=[0,0,0])
    cmd.origin("origin0")
    t_sap = 10  # findMaxDist("tmp_scaffold")/2.0     # approximate half thickness of saposin monomer

    print(f"State of empty/not-empty: {empty}")

    if not empty:
        avXY  = TMdistCheck("tmp_prot", 0.2)
        radXY = avXY / 2.0
        # remove lipids inside pore
        cmd.remove(f"tmp_memb within {radXY} of origin0")
        r_lipHead = 4.7  # Area(POPC) = 65 A^2 => radius = sqrt(A/pi); reasonable estimate
        numLipLayers = 2.0  # number of lipid layers between TM of core and Saposin tmp_scaffold
        inRadius = avXY + numLipLayers*r_lipHead + t_sap  # equatorial position of tmp_scaffold center of mass
    else:
        inRadius = 27.0 # equatorial position of tmp_scaffold center of mass
    outRadius = inRadius + t_sap                          # cut-lipids beyond this distance
    print(f"Inner radius: {inRadius}")
    print(f"Outer radius: {outRadius}")

    n_sym = int(n_sym)
    rotAng = 360. / float(n_sym)

    # chain ID
    id = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    s0 = range(0,n_sym,1)
    # Build symmates with desired rotations
    cmd.rotate("x", "90", "tmp_scaffold")
    for i in s0:
        angle =+ i*rotAng
        cmd.copy(f"seg{i}", "tmp_scaffold")
        cmd.origin("origin0")
        cmd.rotate("y", str(initRotAngle), f"seg{i}")
        cmd.origin("origin0")
        cmd.translate(f"[0,{inRadius},0]", f"seg{i}")
        chn = id[i+len(s0)]
        cmd.alter(f"seg{i}", f"chain = '{chn}'")
        cmd.origin("origin0")
        cmd.rotate("z", str(angle), f"seg{i}")

    # remove lipids beyond border encase by saposins
    cmd.remove(f"br. org and tmp_memb beyond {outRadius} of origin0")
    # remove lipids clashing with tmp_protein core and saposins
    if not empty:
        cmd.remove("br. org and tmp_memb within 0.3 of pol. and not hydro")
    else:
        cmd.remove("br. org and tmp_memb within 0.3 of seg* and not hydro")

    # Combine into a single PyMol object
    if empty == True:
        s = f"{prefixName}empty_{membrane}_{scaffold}_{(int)(initRotAngle)}"
    else:
        s = f"{prefixName}{protein}_{membrane}_{scaffold}_{(int)(initRotAngle)}"
    if refine:
        s = f"{prefixName}{protein}_{membrane}_{scaffold}_{(int)(initRotAngle)}_{(int)(n_sym)}"

    cmd.create(s, f"{protein}, tmp_memb, seg*")
    cmd.save(s + ".pdb", s)

    cmd.delete("tmp_memb")
    cmd.delete("tmp_scaffold")
    cmd.delete("tmp_prot")
    cmd.delete("seg*")
    cmd.delete("origin0")

    return s

def builderNanodisc(protein, membrane, scaffold, prefixName, offset = 0, refine = False):
    """
    builds a MP-nanodisc systems
    scaffold in this case is a double belt of MSP
    """
    # Checking time of builder function execution
    print(f'protein is: {protein}')
    print(f'scaffold is: {scaffold}')
    print(f'membrane is: {membrane}')
    empty = False
    if protein == None: empty = True
    # copies to delete later
    cmd.copy("tmp_scaffold",scaffold) # store initial
    cmd.copy("tmp_memb",membrane) # store initial
    cmd.copy("tmp_prot",protein) # store initial
    center("tmp_memb")
    center("tmp_scaffold")
    cmd.pseudoatom("origin0", pos=[0,0,0])
    cmd.origin("origin0")
    outRadius = findAverDist("tmp_scaffold")
    cmd.translate(f"[0,0,{offset}]", f"tmp_scaffold")
    print(f"Max distance from origin to scaffold in xy plane: {outRadius}")
    # remove lipids beyond border encased by MSP
    cmd.remove(f"org and tmp_memb beyond {outRadius} of origin0")
    print(f"State of empty/not-empty: {empty}")
    # remove lipids clashing with tmp_protein core
    if not empty:
        avXY = TMdistCheck("tmp_prot", 0.2)
        minXY = avXY/2.0
        # remove lipids inside pore
        cmd.remove(f"org and tmp_memb within {minXY} of origin0")
        print(f"Mean distance if TM cross-section in xy plane: {avXY}")

    if empty:
        cmd.remove("org and tmp_memb within 0.4 of tmp_scaffold and not hydro")
        s = f"{prefixName}empty_{membrane}_{scaffold}"
    else:
        cmd.remove("org and tmp_memb within 0.3 of pol. and not hydro")
        s = f"{prefixName}{protein}_{membrane}_{scaffold}"
    if refine: s += str(int(offset))
    cmd.create(s,f"({protein},tmp_scaffold, tmp_memb)")
    cmd.save(s + ".pdb", s)

    cmd.delete("tmp_memb")
    cmd.delete("tmp_scaffold")
    cmd.delete("tmp_prot")
    cmd.delete("origin0")
    return s


def builderDetergent(protein, detergent, prefixName, ang = None, densAng = None, refine = False):
    """
    builds MP - detergent complex using a single detergent molecule
    """

    # Checking object names
    print(f'protein   is: {protein}')
    print(f'detergent is: {detergent}')
    center(detergent)
    #if it is an empty assembly --> build a micelle
    if protein == None:
        # table with common empty micelle parameters
        #https://www.molecularbiophysics.physik.lmu.de/publications/Lipfert_etal_JPCB07.pdf
        if refine:
            radius             =  ang
            numberOfDetergents =  densAng
        else:
            radius             =  findMaxDist(detergent)
            numberOfDetergents = 300
        s = builderMicelle(detergent, radius, numberOfDetergents)
        cmd.save(s + ".pdb", s)
        return s

    cmd.copy("tmp_prot", protein)  # store initial
    cmd.pseudoatom("origin0", pos=[0,0,0])
    # Determine max distance of TM cross-section (xy plane)
    r        = TMdistCheck("tmp_prot", 2.0)
    detR     = findMaxDist(detergent)
    print(f"Max distance if TM cross-section is in a xy plane: {r}")
    print(f"Max distance of detergent : {detR}")
    
    #Shrink detergent along z axis to match ry
    #stretch = r/detR
    #affineStretch(detergent, stretch)
    # find new
    #detR     = detR * stretch
    print(f"Max distance of detergent after shrinking: {detR}")
    # Create a ring of detergents using spherical coordinates
    # TODO: find automatically?
    if refine:
        # stochastic
        #theta = random.sample(range(-ang, ang), 10)
        #phi = random.sample(range(0, 360), densAng)
        # geometrical
        theta = np.arange(-ang, ang + 1, 3)
        phi    = np.arange(0, 361, densAng)
    else:
        # stochastic
        #theta = random.sample(range(-20, 20), 10)
        #phi = random.sample(range(0, 361), 100)
        # geometrical
        theta       = range (-10, 19, 3)
        phi          = range(0, 361, 10)   # find angular step from average density?
    builderCorona(theta, phi, detergent, r, detR)
    # combine components into single PYMOL object
    if prefixName is not "" : s = f"{prefixName}"
    if refine:
        s = f"{protein}_{detergent}_{(int)(ang)}_{(int)(densAng)}"
    else:
        s = f"{protein}_{detergent}"

    #affineStretch("corona", 1.1)
    cmd.create(s, f"{protein}, corona")
    cmd.delete("corona")
    cmd.delete("tmp_prot")
    cmd.delete("origin0")
    cmd.save(s + ".pdb", s)
    return s

def builderMicelle(detergent, r, numberOfDetergents):
    i = 0
    numberOfDetergents = (int)(numberOfDetergents)
    #FIXME: if number of detergents > 360, molecules may clash in space
    if numberOfDetergents > 360:
        x1 = range(-180,180)
        x2 = range(0,360)
        theta = ([random.choice(x1) for _ in range(numberOfDetergents)])
        phi   = ([random.choice(x2) for _ in range(numberOfDetergents)])
    else:
        theta = random.sample(range(-180, 180), numberOfDetergents)
        phi = random.sample(range(0, 360), numberOfDetergents)
    for t, p in zip(theta, phi):
        i += 1
        cmd.copy(f"seg{i}", detergent)
        # randomly sample on a sphere
        cmd.translate(f"[0,{r},0]", f"seg{i}")
        cmd.rotate("x", f"{t}", f"seg{i}")
        cmd.rotate("z", f"{p}", f"seg{i}")
    s = f"micelle_{detergent}_{(int)(r)}_{(int)(numberOfDetergents)}"
    cmd.create(f"{s}", "seg*")
    cmd.delete("seg*")
    #center(f"{s}")
    #cmd.show_as("sticks","org")
    # could be streched if necessary
    # affineStretch(s, 10)
    return s

def builderCorona(theta, fi, detergent, r, detR):
    # Build symmates with desired rotations
    thetaSteps = len(theta)
    angleVer = np.linspace(-90, 90, thetaSteps)
    i = 0
    for t, a in zip(theta, angleVer):
        for f in fi:
            i += 1
            cmd.copy(f"seg{i}", detergent)
            # corona
            cmd.rotate("x", f"{a}", f"seg{i}")
            cmd.translate(f"[0,{r + 0.5*detR*np.cos(np.deg2rad(a))},0]", f"seg{i}")
            cmd.translate(f"[0,0,{(r+detR)*np.sin(np.deg2rad(t))}]", f"seg{i}")
            cmd.rotate("z", f"{f}", f"seg{i}")
    cmd.create("corona", "seg*")
    cmd.delete("seg*")

def builderMembrane(lipid):
    """
    build membrane bilayer from single lipid PDB file
    """
    ###FIXME: add protein to the model and change the modelName if not empty. Otherwise call the model "{protein}_{lipid}"
    cmd.load(lipid+".pdb", "start_lipid")
    cmd.alter("start_lipid", "chain = 'X'")
    cmd.alter("start_lipid", "segi = 'mema'")
    cmd.rotate('x', 90, "start_lipid")
    dmax = findMaxDist("start_lipid")

    # create lipid copies and translate them to new position
    nlip = 20 # number of lipids forming edge of bilayer

    s0 = range(1,nlip,1)
    s1 = range(1,nlip+1,1) # excludes first lipid

    step_x = 0 # translation in x (TODO: automatic determination of spacing without clashes)
    step_y = 7
    step_z = 0
    step_x2 = 7
    step_y2 = 0
    step_z2 = 0

    for i in s1:
        # first column
        cmd.copy(f"lip{i}", "start_lipid") # row of lipids
        cmd.alter(f"lip{i}", f"resi{i}") # change residue numbers
        y = i * step_y
        cmd.translate(f"[{step_x},{y},{step_z}]", f"lip{i}")
        # generate remaining rows/columns in same leaflet
        for j in s0:
            k = int(nlip) * i + j # TODO: general counter to write correct lipid number
            cmd.copy(f"lip{k}", f"lip{i}") # adjacent row of lipids
            cmd.alter(f"lip{k}", f"resi={k}") # change residue numbers
            x2 = j*step_x2
            cmd.translate(f"[{x2},{step_y2},{step_z2}]", f"lip{k}")
        cmd.sort() # sort atom order
    # create second leaflet
    # simple method by creating a single leaflet object:
    cmd.create("mema","(lip*)")
    cmd.delete("lip*")

    cmd.copy("memb","mema")
    cmd.alter("memb", "segi = 'memb'")
    cmd.rotate("x", 180, "memb")
    cmd.translate(f"[0,0,{(-1.0*(dmax + 0.5))}]" , "memb")
    #cmd.color("yellow", "segi = 'mema'")
    #cmd.color("blue", "segi = 'memb'")
    cmd.translate("[3.5,3.5,0]", "memb") # optional shift of single leaflet to avoid aliphatic clashes
    s = f"{lipid}_bilayer"
    cmd.create(s,"(mema,memb)")
    cmd.delete("mema ,memb, start_lipid")
    center(s)
    cmd.save(s + ".pdb", s)
    cmd.reset()
    return s


def builderBicelle(protein, membrane, detergent, prefixName, refine = False, ang = None, densAng = None):
    """
    builds MP - bicelle complex using a membrane and a single detergent molecule
    """
    # Checking object names
    print(f'protein   is: {protein}')
    print(f'membrane  is: {membrane}')
    print(f'detergent is: {detergent}')

    cmd.copy("tmp_prot", protein)  # store initial
    cmd.copy("tmp_memb", membrane)  # store initial
    cmd.copy("tmp_deter", detergent)  # store initial
    center("tmp_memb")
    center("tmp_deter")
    cmd.pseudoatom("origin0", pos=[0, 0, 0])

    center("tmp_prot")
    # Determine max distance of TM cross-section (xy plane)
    r = TMdistCheck("tmp_prot", 2.0)
    detR = findMaxDist("tmp_deter")
    print(f"Max distance if TM cross-section is in a xy plane: {r}")
    print(f"Max distance of detergent : {detR}")

    # Shrink detergent along z axis to match ry
    stretch = r / detR
    affineStretch("tmp_deter", stretch)
    # find new
    detR = detR * stretch
    print(f"Max distance of detergent after shrinking: {detR}")
    # Create a ring of detergents using spherical coordinates
    # FIXME: find automatically?
    if refine:
        # stochastic
        # theta = random.sample(range(-ang, ang), 10)
        # phi = random.sample(range(0, 360), densAng)
        # geometrical
        theta = np.arange(-ang, ang + 1, 5)
        phi = np.arange(0, 360, densAng)
    else:
        # stochastic
        # theta = random.sample(range(-20, 20), 10)
        # phi = random.sample(range(0, 361), 100)
        # geometrical
        theta = range(-20, 21, 5)
        phi = range(0, 361, 5)  # find angular step from average density?

    builderCorona(theta, phi, "tmp_deter", r, detR)
    affineStretch("corona", 1.1)
    # remove lipids inside pore
    cmd.remove(f"tmp_memb within {r/2.0} of origin0")
    cmd.origin("origin0")
    # remove lipids beyond border encased by MSP
    print(f"org and tmp_memb beyond {r} of origin0")
    cmd.remove(f"org and tmp_memb beyond {r} of origin0")
    # remove lipids clashing with tmp_protein core and MSP scaffold and combine into a single PyMol object
    cmd.remove("org and tmp_memb within 0.3 of pol. and not hydro")
    s = f"{prefixName}{protein}_{membrane}_{detergent}"

    cmd.create(s,f"({protein}, corona, tmp_memb)")
    cmd.save(s + ".pdb", s)
    cmd.delete("tmp_prot, tmp_memb, tmp_deter, corona")
    cmd.delete("origin0")
    return s

def TMdistCheck(selection, z):
    """
    determine cross section distances about origin of trans-membrane region of protein PDB
    (distance from origin to protein atoms in xy plane)
    """
    cmd.pseudoatom("origin0", pos=[0,0,10])
    dmax = []
    for ang in range (0, 181, 10):
        dlist = []
        cmd.rotate('z', ang, selection)
        xLine = f"{selection} and z > -{z} and z < {z}  and y > -{z} and y < {z}"
        atoms = cmd.index(xLine)
        if len(atoms) == 0 :
            print(f"No atoms in {xLine}.")
            continue
        # find R in only one cross-section
        for at1 in cmd.index("origin0"):
            for at2 in atoms:
                dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
                dlist.append(dist)
        dmax.append(max(dlist))
    cmd.rotate('z', -180, selection)
    meanXY = np.mean(dmax)
    return meanXY

def refresh():
    cmd.reset()
    cmd.show_as("cartoon", "pol")
    #cmd.show_as("sticks", "org")
    cmd.show_as("spheres", "org")


def findMaxDist(selection):
    """
    finds the longest distance within a molecule
    """
    dlist = []
    for at1 in cmd.index(selection):
        for at2 in cmd.index(selection):
            dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
            dlist.append(dist)
    dmax = max(dlist)
    return dmax

def findAverDist(selection):
    dlist = []
    cmd.pseudoatom("origin0", pos=[0,0,0])
    for at1 in cmd.index("origin0"):
        for at2 in cmd.index(selection)[::30]:
            dist = cmd.get_distance(atom1=at1, atom2=at2, state=0)
            dlist.append(dist)
    return np.mean(dlist)

# parse crysol log file
def parseCrysolLog(logFileName):
    '''Parse Crysol log file, obtain Chi2, Rg and eDens'''
    # will not parse crysol_summary.txt, but the .log file
    # created for each individual run

    chi2 = -9999;
    Rg = -9999;
    eDens = -9999;

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

def writePdb(sel, prefix = ""):
    pdbfn = prefix + sel + ".pdb"
    npdbfn = pdbfn.replace(" or ", "");
    npdbfn = npdbfn.replace(" and ", "");

    try:
        npdbfn = npdbfn.translate(npdbfn.maketrans('', '', " "))
    except ImportError:
        pass
    else:
        try:
            npdbfn = npdbfn.translate(" ")
        except ImportError:
            pass
    print(f"{sel} is saved to {npdbfn}")
    cmd.save(npdbfn, sel)
    return npdbfn

def systemCommand(command, **kwargs):
    status = subprocess.call(command, **kwargs)
    if(0 != status):
        print("WARNING, something went wrong while executing:\n"
                + ' '.join(command))
    return status

def center(selection):
    # center a molecule
    com = cmd.centerofmass(selection, -1.0)
    shift = np.array2string(np.negative(np.array(com)), separator=', ')
    cmd.translate(shift, selection)
    print(f'Shifting {selection} to center')

def affineStretch(selection, stretch):
    # stretch molecule using affine transformations
    stored.altered = []
    cmd.iterate_state(1, selector.process(selection), f"stored.altered.append([x,{stretch}*y,z])")
    cmd.alter_state(1, selection, "(x,y,z) = stored.altered.pop(0)")

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec

def loadModel(filename, pymolObject):
    if (filename):
        try:
            cmd.delete(pymolObject)
            cmd.load(filename+".pdb")
        except Exception as e:
            print(f"For file {filename} and {pymolObject}: {e}")

#Not used at the moment!! for the future
def affine(selection,sh_x=0.1,sh_y=0.1):
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
    return cmd.transform_selection(selection,rotMat)

def curveBilayer(radius=15.0,shift=1.0):
    """ 
    Raise/lower height of central bilayer lipids to create a curved bilayer surface
    """
    shiftTop = shift  # shift in + Z
    shiftBot = -1.0 * float(shift) # shift in -Z 
    cmd.translate([0.0,0.0,shiftTop], "br. segi 'mema' within %s of origin0" %radius) # Adjust height of bilayer
    cmd.translate([0.0,0.0,shiftBot], "br. segi 'memb' within %s of origin0" %radius) # Adjust height of bilayer
    print("Adjusted vertical position of central lipids by %s A"  %shift)
