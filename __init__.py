'''
PyMOL simple membrane builder

The plugin is constructed to facilitate SAXS analysis/modeling of membrane protein systems
(D.Molodenskiy & H.D.T.Mertens 2019-2020 for BioSAXS team)
'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.

import time
from .util import *

from pymol import cmd
# pymol.Qt provides the PyQt5 interface, but may support PyQt4 and/or PySide as well
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi

from sys import platform
viewer = 'primus'
if platform == 'win32' or platform == 'win64':
    viewer = 'primusqt'

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Membrane Builder Plugin', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
mpb = None

def run_plugin_gui():
    '''
    Open plugin dialog
    '''
    global mpb

    if mpb is None:
        mpb = mpbuilder()

    mpb.dialog.show()


class mpbuilder():


    def changeForm(self, assembly):
        if assembly == "detergent":
            self.form.label_6.setEnabled(True)
            self.form.input_filename_prot.setEnabled(True)
            self.form.btn_browse_prot.setEnabled(True)
            self.form.first_angle_label.setEnabled(True)
            self.form.input_rotAng_min.setEnabled(True)
            self.form.input_rotAng_max.setEnabled(True)
            self.form.input_rotAng_step.setEnabled(True)
            self.form.first_angle_label.setText(r'<html><head/><body><p><span style=" font-size:10pt;">Max Polar Angle </span><span style=" font-size:10pt; font-weight:600;">ω</span><span style=" font-size:10pt;"> (deg.)</span></p></body></html>')
            self.form.second_angle_label.setText(r'<html><head/><body><p><span style=" font-size:10pt;">Density Angle </span><span style=" font-size:10pt; font-weight:600;">Δφ</span><span style=" font-size:10pt;"> (deg.)</span></p></body></html>')
            self.form.label_7.setText("Detergent")
            self.form.btn_browse_scaffold.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.filename_scaffold_label.setEnabled(False)
            self.form.copies_scaffold_label.setEnabled(False)
            self.form.input_copies.setEnabled(False)
            self.form.input_filename_scaffold.setEnabled(False)
            self.form.btn_browse_scaffold.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.input_scaffold_number_min.setEnabled(False)
            self.form.input_scaffold_number_max.setEnabled(False)
            self.form.input_scaffold_number_step.setEnabled(False)
            self.form.second_angle_label.setEnabled(True)
            self.form.input_rotAng_min_2.setEnabled(True)
            self.form.input_rotAng_max_2.setEnabled(True)
            self.form.input_rotAng_step_2.setEnabled(True)
            self.form.label_2.setEnabled(False)
            self.form.input_rotAng.setEnabled(False)
            self.form.checkBox_prebuild_bilayer.setEnabled(False)
            self.form.checkBox_emptyAssembly.setEnabled(True)
        elif assembly == "salipro":
            self.form.label_6.setEnabled(True)
            self.form.input_filename_prot.setEnabled(True)
            self.form.btn_browse_prot.setEnabled(True)
            self.form.first_angle_label.setEnabled(True)
            self.form.input_rotAng_min.setEnabled(True)
            self.form.input_rotAng_max.setEnabled(True)
            self.form.input_rotAng_step.setEnabled(True)
            self.form.first_angle_label.setText("Scaffold Rotation Angle (deg.)")
            self.form.label_7.setText("Lipid Bilayer")
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.input_copies.setEnabled(True)
            self.form.number_scaffold_label.setEnabled(True)
            self.form.filename_scaffold_label.setEnabled(True)
            self.form.copies_scaffold_label.setEnabled(True)
            self.form.filename_scaffold_label.setText("Scaffold")
            self.form.input_filename_scaffold.setEnabled(True)
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.number_scaffold_label.setEnabled(True)
            self.form.input_scaffold_number_min.setEnabled(True)
            self.form.input_scaffold_number_max.setEnabled(True)
            self.form.input_scaffold_number_step.setEnabled(True)
            self.form.second_angle_label.setEnabled(False)
            self.form.input_rotAng_step_2.setEnabled(False)
            self.form.input_rotAng_min_2.setEnabled(False)
            self.form.input_rotAng_max_2.setEnabled(False)
            self.form.label_2.setEnabled(True)
            self.form.input_rotAng.setEnabled(True)
            self.form.checkBox_prebuild_bilayer.setEnabled(True)
            self.form.checkBox_emptyAssembly.setEnabled(True)
        elif assembly == "nanodisc":
            self.form.label_6.setEnabled(True)
            self.form.input_filename_prot.setEnabled(True)
            self.form.btn_browse_prot.setEnabled(True)
            self.form.first_angle_label.setEnabled(True)
            self.form.first_angle_label.setText("Offset along X-axis (A)")
            self.form.input_rotAng_min.setEnabled(True)
            self.form.input_rotAng_max.setEnabled(True)
            self.form.input_rotAng_step.setEnabled(True)
            self.form.label_7.setText("Lipid Bilayer")
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.input_copies.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.filename_scaffold_label.setEnabled(True)
            self.form.copies_scaffold_label.setEnabled(False)
            self.form.filename_scaffold_label.setText("Scaffold")
            self.form.input_filename_scaffold.setEnabled(True)
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.input_scaffold_number_min.setEnabled(False)
            self.form.input_scaffold_number_max.setEnabled(False)
            self.form.input_scaffold_number_step.setEnabled(False)
            self.form.second_angle_label.setText("Offset along Y-axis (A)")
            self.form.second_angle_label.setEnabled(True)
            self.form.input_rotAng_step_2.setEnabled(True)
            self.form.input_rotAng_min_2.setEnabled(True)
            self.form.input_rotAng_max_2.setEnabled(True)
            self.form.input_rotAng.setEnabled(False)
            self.form.checkBox_prebuild_bilayer.setEnabled(True)
            self.form.checkBox_emptyAssembly.setEnabled(True)
        elif assembly == "bilayer":
            self.form.label_6.setEnabled(False)
            self.form.input_filename_prot.setEnabled(False)
            self.form.btn_browse_prot.setEnabled(False)
            self.form.first_angle_label.setEnabled(False)
            self.form.input_rotAng_min.setEnabled(False)
            self.form.input_rotAng_max.setEnabled(False)
            self.form.input_rotAng_step.setEnabled(False)
            self.form.label_7.setText("A single molecule")
            self.form.btn_browse_scaffold.setEnabled(False)
            self.form.input_copies.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.filename_scaffold_label.setEnabled(False)
            self.form.copies_scaffold_label.setEnabled(False)
            self.form.input_filename_scaffold.setEnabled(False)
            self.form.btn_browse_scaffold.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.input_scaffold_number_min.setEnabled(False)
            self.form.input_scaffold_number_max.setEnabled(False)
            self.form.input_scaffold_number_step.setEnabled(False)
            self.form.second_angle_label.setEnabled(False)
            self.form.input_rotAng_step_2.setEnabled(False)
            self.form.input_rotAng_min_2.setEnabled(False)
            self.form.input_rotAng_max_2.setEnabled(False)
            self.form.label_2.setEnabled(False)
            self.form.input_rotAng.setEnabled(False)
            self.form.checkBox_prebuild_bilayer.setEnabled(False)
            self.form.checkBox_emptyAssembly.setEnabled(False)
        elif assembly == "bicelle":
            self.form.label_6.setEnabled(True)
            self.form.input_filename_prot.setEnabled(True)
            self.form.btn_browse_prot.setEnabled(True)
            self.form.first_angle_label.setEnabled(False)
            self.form.input_rotAng_min.setEnabled(False)
            self.form.input_rotAng_max.setEnabled(False)
            self.form.input_rotAng_step.setEnabled(False)
            self.form.label_7.setText("Lipid Bilayer")
            self.form.input_copies.setEnabled(False)
            self.form.number_scaffold_label.setEnabled(True)
            self.form.filename_scaffold_label.setText("Detergent")
            self.form.filename_scaffold_label.setEnabled(True)
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.copies_scaffold_label.setEnabled(False)
            self.form.input_filename_scaffold.setEnabled(True)
            self.form.btn_browse_scaffold.setEnabled(True)
            self.form.number_scaffold_label.setEnabled(False)
            self.form.input_scaffold_number_min.setEnabled(False)
            self.form.input_scaffold_number_max.setEnabled(False)
            self.form.input_scaffold_number_step.setEnabled(False)
            self.form.second_angle_label.setEnabled(False)
            self.form.input_rotAng_step_2.setEnabled(False)
            self.form.input_rotAng_min_2.setEnabled(False)
            self.form.input_rotAng_max_2.setEnabled(False)
            self.form.label_2.setEnabled(False)
            self.form.input_rotAng.setEnabled(False)
            self.form.checkBox_prebuild_bilayer.setEnabled(True)
            self.form.checkBox_emptyAssembly.setEnabled(False)

    def change_assembly(self):
        # Change active buttons
        # type of protein-membrane assembly
        assemblyType = self.form.input_type.currentText()
        print('Assembly type is changed to: {}'.format(assemblyType))
        self.changeForm(assemblyType)

    def preOriProt(self):
        '''
        User defined orientation of protein true/false
        '''
        if self.form.checkBox_3.isChecked():
            #print(f"Using pre-oriented membrane protein {self.protName} for model building ...")
            print("Using pre-oriented membrane protein {} for model building ...".format(self.protName))
        elif self.protName == None:
            print("No protein to pre-align")
        else:
            #print(f"Centering and aligning membrane protein {self.protName} for model building ...")
            print("Centering and aligning membrane protein {} for model building ...".format(self.protName))
            center(self.protName)

    def emptyAssembly(self):
        assemblyType = self.form.input_type.currentText()
        if self.form.checkBox_emptyAssembly.isChecked():
          print("No transmembrane protein, empty assembly will be generated ...")
          self.form.label_6.setEnabled(False)
          self.form.input_filename_prot.setEnabled(False)
          self.form.btn_browse_prot.setEnabled(False)
          self.protName = None
          if assemblyType == "detergent":
            self.form.first_angle_label.setText("Radius-vector (A)")
            self.form.second_angle_label.setText("Number of detergent molecules")
        else:
          print("Transmembrane protein {} will be used for for model building ...".format(self.protName))
          if assemblyType == "detergent":
              self.form.first_angle_label.setText("Maximum Polar Angle Theta  (deg.)")
              self.form.second_angle_label.setText("Density Angle Phi (deg.)")

    def preBuildBilayer(self):
        '''
        activate the button and store the variable
        '''
        if self.form.checkBox_prebuild_bilayer.isChecked():
            print("Using pre-built membrane {}".format(self.membName))
            self.buildMemb = False
        else:
            print("Using a single molecule to build a membrane")
            self.buildMemb = True

    # CRYSOL predict scattering on single assembly
    def run_crysol_predict(self):
        """Prediction of model scattering"""
        #def predcrysol(crycalc, models, prefix="tmp", param=" "):
        if (self.modelName):
            df = predcrysol(self.modelName, "yes")
            if os.path.exists(df):
                print('Theoretical SAXS profile generated')
                systemCommand([viewer, df])
                self.prediction = df
        else:
            print("pdb file is missing. Model Name: {}".format(self.modelName))

    # CRYSOL fit run on single assembly
    def run_crysol_fit(self):
        """Calculation of model fit"""
        if (self.modelName != "" and self.dataName != ""):
            fitcrysol(self.modelName,self.dataName, "yes", True)
        else:
            print('pdb or SAXS data file is missing!')
            print('model file: {} data file: {}'.format(self.modelName, self.dataName))

    def run_refinement(self):
        """Finds the model that fits SAXS data"""
        if self.membName == None:
            print("Please provide a PDB file!")
            return
        if self.dataName == None:
            print("Please provide SAXS dat file!")
            return
        seconds_init = time.time()
        assemblyType = self.form.input_type.currentText()
        print('Assembly is: {}'.format(assemblyType))
        prefixName = self.form.output_filename_prefix.text()
        # center protein
        self.preOriProt()

        rot_min_ang = self.form.input_rotAng_min.value()
        rot_max_ang = self.form.input_rotAng_max.value() + 1
        rot_step_ang = self.form.input_rotAng_step.value()
        # type of protein-membrane assembly
        if assemblyType == "detergent":
            dens_min_ang = self.form.input_rotAng_min_2.value()
            dens_max_ang = self.form.input_rotAng_max_2.value() + 1
            dens_step_ang = self.form.input_rotAng_step_2.value()
            bestModel, fit = crysolRefinementDetergent(rot_min_ang, rot_max_ang,rot_step_ang, \
                                                       dens_min_ang,dens_max_ang,dens_step_ang, \
                                                       self.protName, self.membName, self.dataName,\
                                                       prefixName)


        elif assemblyType == "salipro":
            cmd.reset()
            scaffold_min  = self.form.input_scaffold_number_min.value()
            scaffold_max  = self.form.input_scaffold_number_max.value() + 1
            scaffold_step = self.form.input_scaffold_number_step.value()
            if (self.buildMemb):
                self.membName = builderMembrane(self.membName)
                self.form.input_filename_lip.setText(self.membName)
            refresh()
            bestModel, fit = crysolRefinementSalipro(rot_min_ang, rot_max_ang,rot_step_ang, \
                                                     scaffold_min,scaffold_max,scaffold_step, \
                                                     self.protName, self.membName, self.scafName,self.dataName,\
                                                     prefixName)
        elif assemblyType == "nanodisc":
            cmd.reset()
            y_min  = self.form.input_rotAng_min.value()
            y_max  = self.form.input_rotAng_max.value() + 1
            y_step = self.form.input_rotAng_step.value()
            if (self.buildMemb):
                self.membName = builderMembrane(self.membName)
                self.form.input_filename_lip.setText(self.membName)
            refresh()
            bestModel, fit = crysolRefinementNanodisc(rot_min_ang, rot_max_ang,rot_step_ang, \
                                                      y_min, y_max, y_step, \
                                                      self.protName, self.membName, self.scafName, self.dataName,\
                                                      prefixName)
        else:
            print("Refinement is not supported for assembly type {}".format(assemblyType))
            return
        self.modelName = bestModel
        self.fit = fit
        seconds_tmp = time.time()
        refresh()
        t = int((seconds_tmp - seconds_init))
        print("{:d} seconds consumed.".format(t))
        systemCommand([viewer, self.fit])

    def run_build(self):
        #callback for the "Build" button
        seconds_init = time.time()
        # get form data
        prefixName = self.form.output_filename_prefix.text()
        assemblyType = self.form.input_type.currentText()
        print('Assembly is: {}'.format(assemblyType))
        print('Building model...')
        # center protein if needed
        self.preOriProt()
        #delete old model
        cmd.delete(self.modelName)
        # executions depends on the assembly type
        if assemblyType == "detergent":
            # execute detergent builder
            self.modelName = builderDetergent(self.protName, self.membName, prefixName)
        elif assemblyType == "salipro":
            #  execute salipro builder
            rotAng = self.form.input_rotAng.value()
            numScaffoldCopies = self.form.input_copies.value()
            print('Number of scaffold copies: {}'.format(numScaffoldCopies))
            if (self.buildMemb):
                self.membName = builderMembrane(self.membName)
                self.form.input_filename_lip.setText(self.membName)
            self.modelName = builderSalipro(self.protName, self.scafName, self.membName, prefixName,
                                              numScaffoldCopies, rotAng)

        elif assemblyType == "nanodisc":
            # execute nanodisc builder
            # build bilayer if check box is activated
            if (self.buildMemb):
                self.membName = builderMembrane(self.membName)
                self.form.input_filename_lip.setText(self.membName)
            #  execute builder
            self.modelName = builderNanodisc(self.protName, self.membName, self.scafName, prefixName)

        elif assemblyType == "bilayer":
            # execute detergent builder
            self.membName = builderMembrane(self.membName)
            self.form.input_filename_lip.setText(self.membName)

        elif assemblyType == "bicelle":
            # execute bicelle builder
            # build bilayer if check box is activated
            if (self.buildMemb):
                self.membName = builderMembrane(self.membName)
                self.form.input_filename_lip.setText(self.membName)
            #  execute builder
            self.modelName = builderBicelle(self.protName, self.membName, self.scafName, prefixName)

        refresh()
        print("Model name is {}".format(self.modelName))
        seconds_tmp = time.time()
        t = int((seconds_tmp - seconds_init))
        print("{:d} seconds consumed.".format(t))

    def browse_filename_data(self, _str):
        filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open file', os.getcwd(), "dat files (*.dat);;")[0]
        if filename:
            f = os.path.relpath(filename)
            self.form.input_filename_data.setText(f)
            self.dataName = f

    def browse_filename_prot(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open protein file', os.getcwd(),
                                                         "PDB files (*.pdb *.cif *.ent);;")[0]
        filename = filename[:-4]  # remove .pdb
        loadModel(filename, self.protName)
        f = os.path.basename(filename)
        self.form.input_filename_prot.setText(f)
        self.protName = f
        refresh()

    def browse_filename_lip(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open lipid/detergent file', os.getcwd(),
                                                         "PDB files (*.pdb *.cif *.ent);;")[0]
        filename = filename[:-4]  # remove .pdb
        loadModel(filename, self.membName)
        f = os.path.basename(filename)
        self.form.input_filename_lip.setText(f)
        self.membName = f
        refresh()

    def browse_filename_scaffold(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(None, 'Open scaffold file', os.getcwd(),
                                                         "PDB files (*.pdb *.cif *.ent);;")[0]
        filename = filename[:-4] #remove .pdb
        loadModel(filename, self.scafName)
        f = os.path.basename(filename)
        self.form.input_filename_scaffold.setText(f)
        self.scafName = f
        refresh()

    def __init__(self):
        # create a new Window
        self.dialog = QtWidgets.QDialog()
        # populate the Window from our *.ui file which was created with the Qt Designer
        uifile = os.path.join(os.path.dirname(__file__), 'widgetMPbuilder.ui')
        self.form = loadUi(uifile, self.dialog)

        # open always on detergents
        self.changeForm("detergent")
        self.protName = self.membName = self.scafName = self.dataName = None
        # name of the model
        self.modelName    = ""
        # prefix
        self.prefixName = ""
        # crysol int prediction
        self.prediction = ""
        #crysol fit
        self.fit        = ""
        #result of the fit
        self.fitResult = {'chi2': -9999, 'Rg': -9999, 'eDens': -9999}
        #whether to build a membrane or use prebuilt one
        self.buildMemb = False
        # hook up button callbacks
        self.form.button_build.clicked.connect(self.run_build)
        self.form.button_refinement.clicked.connect(self.run_refinement)
        self.form.button_crysol_fit.clicked.connect(self.run_crysol_fit)  # single fit of assembly to data
        self.form.button_crysol_predict.clicked.connect(self.run_crysol_predict)  # predict theory
        #self.form.button_crysol_predict.clicked.connect(lambda: self.run_crysol_predict(modelName))

        self.form.btn_browse_prot.clicked.connect(self.browse_filename_prot)
        self.form.btn_browse_lip.clicked.connect(self.browse_filename_lip)
        self.form.btn_browse_scaffold.clicked.connect(self.browse_filename_scaffold)
        self.form.btn_browse_data.clicked.connect(self.browse_filename_data)

        self.form.button_close.clicked.connect(self.dialog.close)
        
        # change GUI depending upon the choice of assembly type
        self.form.input_type.currentIndexChanged.connect(self.change_assembly)
        # emptyAssembly = assemble complex without transmembrane protein component
        self.form.checkBox_emptyAssembly.toggled.connect(self.emptyAssembly)
        # preOriProt = pre-oriented protein
        self.form.checkBox_3.toggled.connect(self.preOriProt)
        self.form.checkBox_prebuild_bilayer.toggled.connect(self.preBuildBilayer)
