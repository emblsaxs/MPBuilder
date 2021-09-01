# MPBuilder
MPBuilder provides building capability for protein-detergent, bicelle, and lipid-scaffold (saposin nanoparticles, nanodiscs) complexes and links this to the ATSAS software package modules for model refinement and validation against the SAXS data.
Below is the example GUI snapshot for detergent complexes.

![alt text](https://github.com/emblsaxs/MPBuilder/blob/main/gui.png?raw=true)

The membrane protein assemblies available for modeling by MPBuilder are presented below.

![alt text](https://github.com/emblsaxs/MPBuilder/blob/main/fig1mm.png?raw=true)

## Installation ##
* Install an Open-Source PyMOL version (see below).
* Make sure ATSAS is installed and working.
* Download the *mpbuilder* directory from https://github.com/emblsaxs/MPBuilder. If ATSAS 3.0.4 or later version is used the MPBuilder folder is automatically copied from github at "C:\Program Files\ATSAS-3.xxx.yyy\bin\python\site-packages\pymolplugins\mpbuilder" upon installation.
* Make a zip archive from the downloaded directory, e.g. mpbuilder.zip
* Create or modify $HOME/.pymolrc.pml by adding the following line:
  > os.environ["PATH"] += os.pathsep + "/xxx/yyy/bin:"
  
  where /xxx/yyy is the path to your local ATSAS installation.
* Start PyMOL
* Go to _Plugin_->_Plugin Manager_->_Install New Plugin_->_Install from local file_
* Browse to the *mpbuilder.zip*, select it, and click _Open_

For more details:  
  * https://pymolwiki.org/index.php/Plugins
  * https://pymolwiki.org/index.php/Pymolrc

## Please use PyMOL 2.0 or the Open-Source version ##

### PyMOL 2.0 ###
* PyMOL version 2.0 has been tested and works properly with MPBuilder.
* It can be downloaded and installed free of charge from:
	* https://pymol.org/2/#download
* It is available for Linux, Mac and Windows.


### Open-Source PyMOL ###

Open-source PyMOL installations:

* Linux - on your terminal type: > sudo apt-get install pymol
  * apt-get is for Debian/Ubuntu, please adapt depending on your distro, more details at:
  * https://pymolwiki.org/index.php/Linux_Install#Open-Source_PyMOL_in_Linux_Distros
  
* Mac OS X  - installation command: > brew install homebrew/science/pymol
  * For more details: https://pymolwiki.org/index.php/MAC_Install#Pre-compiled
  * (please avoid hybrid and other versions, these may not work properly)
  
* Windows - please follow instructions at:
  * https://pymolwiki.org/index.php/Windows_Install#Open-Source_PyMOL
  
  
In case of doubts, bugs, problems or comments please write to:
atsas@embl-hamburg.de
