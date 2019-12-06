..  -*- coding: utf-8 -*-

Getting Started with Scopy
==========================
This document is intended to provide an overview of how one can use the Scopy functionality from Python. If you find mistakes, or have suggestions for improvements, please either fix them yourselves in the source document (the .py file) or send them to the mailing list: oriental-cds@163.com and yzjkid9@gmail.com.


Installing the Scopy package
-----------------------------
PyBioMed has been successfully tested on Linux and Windows systems under python3 enviroment.

Dependencies
~~~~~~~~~~~~
.. code-block:: python

	RDkit>=2019.03.1
	Numpy

Install from source
~~~~~~~~~~~~~~~~~~~
>>> git clone git@github.com:kotori-y/Scopy.git && cd scopy
>>> [sudo] python setup.py install

PyPI
~~~~
>>> 

Conda
~~~~~
>>> 

Drug-likeness Filter
--------------------
The :mod:`druglikeness` provides the method to analyse the physicochemical (PC) properties and filter compounds based on PC-derived rules. 

Calculating PC properties
~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`druglikeness.molproperty` module provides the tool to calculate PC properties.

>>> import os
>>> from rdkit import Chem
>>> from scopy import ScoConfig
>>> from scopy.druglikeness import molproperty

>>> mol = Chem.MolFromMolFile(os.path.join(ScoConfig.DemoDir,'mol.sdf'))
>>> mol
<rdkit.Chem.rdchem.Mol object at 0x0000020879B4E120>

You could calculate different properties respectively

>>> MW = molproperty.CalculateMolWeight(mol)
>>> MW
229.24
>>> logP = molproperty.CalculateLogP(mol)
>>> logP
2.35
>>> nHD = molproperty.CalculateNumHDonors(mol)
>>> nHD
2

You could also calculate total(43) properties at once time

>>> props = molproperty.GetProperties(mol)
>>> props
{'MW': 229.24, 'Vol': 235.2, 'Dense': 0.97, 'fChar': 0, 'nBond': 18, 'nAtom': 28, 'nHD': 2, 'nHA': 4, 'nHB': 6, 'nHet': 4, 'nStero': 0, 'nHev': 17, 'nRot': 1, 'nRig': 14, 'nRing': 2, 'logP': 2.35, 'logD': 0.8670776309515202, 'pKa': -5.931602224875785, 'logSw': -2.95, 'ab': 'base', 'MR': 64.8, 'tPSA': 69.89, 'AP': 0.35, 'HetRatio': 0.31, 'Fsp3': 0.08, 'MaxRing': 6, 'QEDmean': 0.73, 'QEDmax': 0.7, 'QEDnone': 0.79, 'SAscore': 2.96, 'NPscore': 0.49, 'nSingle': 8, 'nDouble': 4, 'nTriple': 0, 'nC': 13, 'nB': 0, 'nF': 0, 'nCl': 0, 'nBr': 0, 'nI': 0, 'nP': 0, 'nS': 0, 'nO': 3, 'nN': 1}

When calculating the property of multiple molecules, in addition to repeatedly calling the function in :mod:`druglikeness.molproperty`, you can also use :mod:`druglikeness.druglikeness` module, which is more time-saveing since using multiprocessing.

>>> from scopy.druglikeness import druglikeness
>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir, '760.sdf'))
>>> mols = (mol for mol in suppl if mol)
>>> pc = druglikeness.PC_properties(mols=mols, n_jobs=4)
>>> res = pc.CalculateMolWeight()
>>> len(res)
760
>>> type(res)
<class 'list'>
>>> res[:10]
[256.26, 288.25, 182.17, 578.53, 592.55, 286.24, 432.38, 270.24, 448.38, 578.52]

The function return a `list`. Parameter `mols` should be an iterable object (i.g. `list`, `tuple` or `generator`) and `n_jobs` is the number of CPUs to use to do the computation, -1 means using all processors.

Filtering molecule under PC-derived rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`druglikeness.rulesfilter` module provides the tool to analyse PC properties

>>> from scopy.druglikeness import rulesfilter
>>> res = rulesfilter.CheckLipinskiRule(mol)
>>> res
{'Disposed': 'Accepted', 'nViolate': 0}

The function return a `dict`, the field :mod:`Disposed` represents compound state after filter applied (**Rejected** meant the compound rejected by filter, **Accepted** for accepted); :mod:`nViolate` represents the number of PC property violated by compound.

In above example, the compound do not violate any property limited in Lipinski Rule and its status is 'Accepted'.

Besides, the specific value of each propety would be returned if the :mod:`detail` has been set as :mod:`True`.

>>> res = rulesfilter.CheckLipinskiRule(mol, detail=True)
>>> res
{'MW': 229.24, 'logP': 2.35, 'nHD': 2, 'nHA': 4, 'Disposed': 'Accepted', 'nViolate': 0}

You also could customize the filter by your experience

>>> prop_kws = {'MW':[100,500], 'nHB':[5,10], 'QEDmean':[0.8,None]}
>>> res = rulesfilter.Check_CustomizeRule(mol, prop_kws=prop_kws, detail=True)
>>> res
{'MW': 229.24, 'nHB': 6, 'QEDmean': 0.73, 'nViolate': 1, 'VioProp': ['QEDmean']}

The customize rule should be a `dict`, key of `dict` is abbreviation name of property and value is the limited range.

Samely, :mod:`druglikeness.druglikeness` could also be used to analyse multiple molecules, instead of repeatly calling function in `druglikeness.rulesfilter`, to save time

>>> rule = druglikeness.PC_rules(mols,n_jobs=4,detail=True)
>>> res = rule.CheckLipinskiRule()
>>> len(res)
760
>>> type(res)
<class 'list'>
>>> res[:3]
[{'MW': 256.26, 'logP': 2.83, 'nHD': 3, 'nHA': 3, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 288.25, 'logP': 2.79, 'nHD': 5, 'nHA': 5, 'Disposed': 'Accepted', 'nViolate': 0},
 {'MW': 182.17, 'logP': -3.59, 'nHD': 6, 'nHA': 6, 'Disposed': 'Accepted', 'nViolate': 1}]

Structure Alert Filter
----------------------
The :mod:`structure_alert` module provides the tool to search for the presence of toxicophores and flag unwanted reactive chemical groups, where both toxicophores and unwanted reactive chemical groups are encoded by SMARTS.

Screening a single molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`FilterWithSmarts` module provides the tool to screen a molecule. 

>>> from scopy.structure_alert import FilterWithSmarts
>>> mol = Chem.MolFromMolFile(os.path.join(ScoConfig.DemoDir,'PAINS.sdf'))
>>> mol
<rdkit.Chem.rdchem.Mol object at 0x000001E75CE60580>

In here, the PAINS Filter be used to screen a molecule

>>> res = FilterWithSmarts.Check_PAINS(mol)
>>> res
{'Disposed': 'Rejected', 'Endpoint': 'Pains'}

The function return a `dict`, the field :mod:`Disposed` represents compound state after filter applied (**Rejected** meant the compound rejected by filter, **Accepted** for accepted); :mod:`Endpoint` represents the which filter to be used.

Besides, the more specific information would be returned, if the :mod:`detail` has been set as :mod:`True` (defaults to :mod:`False`)

>>> res = FilterWithSmarts.Check_PAINS(mol, detail=True)
>>> res
{'Disposed': 'Rejected', 'MatchedAtoms': [((3, 2, 1, 0, 15, 16, 13, 12),)], 
 'MatchedNames': ['Quinone_A'], 'Endpoint': 'Pains'}

The result reveals the compound rejected by PAINS Filter, since the compound has the substructure named 'Quinone_A' which contained in PAINS Filter, more further, the No.3, No.2, No.1, No.0, No.15, No.16, No.13 and No.12 atom constructing this substructure.

Screening multi-molecule
~~~~~~~~~~~~~~~~~~~~~~~~
In reality, we trend to screen the compund library rather than sinlgle molecule. The :mod:`SmartsFilter` module provides the tool to screen multi-molecule 

>>> from scopy.structure_alert import SmartsFilter
>>> import time

>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir,'760.sdf'))
>>> mols = (mol for mol in suppl if mol)

>>> start = time.process_time()
>>> F = SmartsFilter.Filter(mols, n_jobs=4, detail=True)
>>> res = F.Check_PAINS()
>>> end = time.process_time()
>>> print('Time cost: {}\nSamples: {}'.format((end-start), len(res)))
Time cost: 1.578125
Samples: 760

In above example, the PAINS Filter used to screen a library which contains 760 molecules.
Under using four theardings, the process has cost about 1.6s.

The function return a `list`

>>> res[0]
{'Disposed': 'Accepted', 'MatchedAtoms': ['-'], 'MatchedNames': ['-'], 'Endpoint': 'Pains'}
>>> res[207]
{'Disposed': 'Rejected', 'MatchedAtoms': [((7, 16, 15, 17, 18, 19, 20, 21, 14),)], 'MatchedNames': ['Mannich_A'], 'Endpoint': 'Pains'}

Fingerprint Calculator
----------------------
The :mod:`fingerprint` module provides the tool to compute fingerprints retrieved from fragments.

Comupting fingerprint of single molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There several modules in :mod:`fingerprint` provided to compute corresponding fingerprint 

>>> from scopy.fingerprint import maccs
>>> mol = Chem.MolFromMolFile(os.path.join(ScoConfig.DemoDir,'sid2244.sdf'))
>>> ma = maccs.MACCS()
>>> fp = ma.CalculateMACCS(mol)
>>> fp
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0,
1, 1, 0, 1, 1, 1, 1, 0]

>>> from scopy.fingerprint import efg
>>> efg_fp = efg.EFG(useCount=False)
>>> fp = efg_fp.CalculateEFG(mol)
>>> len(fp)
853

The more details of EFG fingerprint in `Salmina, Elena, Norbert Haider and Igor Tetko (2016)`_

There are 7 types of fingerprint are implemnted in Scopy: MACCS, EFG, EState, Morgan, GhoseCrippen, Daylight and PubChem.

.. _Salmina, Elena, Norbert Haider and Igor Tetko (2016):
	https://www.mdpi.com/1420-3049/21/1/1

Comupting fingerprint of multi-molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Samely, :mod:`fingerptints` can compute multi-molecule

>>> suppl = Chem.SDMolSupplier(os.path.join(ScoConfig.DemoDir,'760.sdf'))
>>> mols = (mol for mol in suppl if mol)
>>> fps = fingerprints.CalculateECFP(mols, radius=2, nBits=1024, n_jobs=4)
>>> fps
array([[0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 1, 0, 0],
       [0, 1, 0, ..., 0, 0, 0],
       ...,
       [0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0],
       [0, 0, 0, ..., 0, 0, 0]])
>>> fps.shape
(760, 1024)

The function return a `numpy.ndarray`. In above example, 760 molecules are calculted ECFP4 fingerprints with 1024 bits.

Screening Visualizer
--------------------

PC properties visualizer
~~~~~~~~~~~~~~~~~~~~~~~~

PC-drivered rules visualizer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fragment visualizer
~~~~~~~~~~~~~~~~~~~

Fingerprint visualizer
~~~~~~~~~~~~~~~~~~~~~

Molecular Pretreater
-------------------





.. figure:: /image/CPI.png
	:width: 400px
	:align: center
	
	The calculation process for chemical-protein interaction descriptors.



























