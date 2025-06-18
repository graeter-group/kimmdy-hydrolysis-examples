from pathlib import Path
from src.utils import gather_0_dists_and_forces

root = Path('/hits/fast/mbm/buhrjk/phd/kimmdy-examples').resolve()
example = root / "examples/collagen-hydrolysis"
systems = ["collagen"]
system = systems[0]

gather_0_dists_and_forces(example, system, dt=0)
