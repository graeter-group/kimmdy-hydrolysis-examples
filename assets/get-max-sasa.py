import MDAnalysis as mda
from mdakit_sasa.analysis.sasaanalysis import SASAAnalysis

u = mda.Universe('./assets/gly.pdb')

ix_cc = 8
ix_n = 13
sasa = SASAAnalysis(u, select=f"index {ix_cc}")
sasa.run(verbose=False)
print(sasa.results.total_area)
