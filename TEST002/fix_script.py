import
os
import
shutil
with open('app/models/adsorption_calculator_fixed.py', 'w', encoding='utf-8') as f:
    f.write('import os\\nimport numpy as np\\nfrom pymatgen.core.structure import Structure, Molecule\\nfrom pymatgen.core.sites import PeriodicSite\\nfrom pymatgen.analysis.adsorption import AdsorbateSiteFinder\\nfrom pymatgen.io.cif import CifWriter\\nfrom pymatgen.symmetry.analyzer import SpacegroupAnalyzer\\nfrom pymatgen.core.surface import SlabGenerator')
