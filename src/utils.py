import os
import random
import subprocess as sp
from contextlib import contextmanager
from pathlib import Path
from string import Template

import pandas as pd
from kimmdy.plugin_utils import bondstats_to_csv, calculate_bondstats
from kimmdy.topology.topology import Topology
from plotnine import ggplot

from src.constants import CLUSTER, STRIDE, T
from src.style import DPI, double_column


def show(p: ggplot, s: float = 1.5):
    path = "tmp/tmp.png"
    os.remove(path) if os.path.exists(path) else None
    p.save(
        path,
        width=double_column,
        height=double_column,
        dpi=DPI,
        verbose=False,
    )


@contextmanager
def pushd(path):
    """Context manager to change directory and return to previous directory."""
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def sh(cmd):
    """Run shell command."""
    sp.run(cmd, text=True, shell=True, check=True)


def read_xvg(
    filename: str | Path,
    columns: list[str] | None = None,
    fields: list[int] | None = None,
) -> pd.DataFrame:
    ls = []
    with open(filename, "r") as f:
        if columns is None:
            # extract column names from the file
            columns = []
            l = f.readline()
            while l.startswith("#") or l.startswith("@"):
                # x column
                if "xaxis" in l:
                    cname = l.split('"')[1]
                    columns.append(cname)
                # y columns
                if l.startswith("@ s"):
                    cname = l.split('"')[1]
                    columns.append(cname)

                l = f.readline()

        for l in f:
            if l.startswith("@") or l.startswith("#"):
                continue
            l = l.split()
            if fields is not None:
                l = [l[i] for i in fields]
            # silently ignore incomplete lines
            # usually the last line
            if len(l) == len(columns):
                ls.append(l)
    return pd.DataFrame(ls, columns=columns, dtype=float)  # pyright: ignore


def write_topology_info(example, system, stride=STRIDE):
    """Write information about backbone bonds and generate a plumed.dat input files."""
    top = Topology.from_path(f"{example}/assets/{system}.top")

    # backbone bonds
    backbone_atom_names = ["C", "N", "CA"]

    # just the peptide bonds (to match up with keys for hydrolysis recipe steps)
    # see get_peptide_bonds_from_top
    ls_peptide = ["key,ai,aj,c,n"]

    ls_backbone = ["key,ai,aj,c,n"]
    backbone = []
    i = 0
    i_peptide = 0
    for bond in top.bonds.keys():
        a = top.atoms[bond[0]]
        b = top.atoms[bond[1]]
        if a.residue in ["NME", "ACE"] or b.residue in ["NME", "ACE"]:
            continue
        if a.atom == "C" and b.atom == "N":
            ls_peptide.append(
                f"{i_peptide},{bond[0]},{bond[1]},{a.residue},{b.residue}"
            )
            i_peptide += 1
        if a.atom in backbone_atom_names and b.atom in backbone_atom_names:
            backbone.append((bond[0], bond[1]))
            ls_backbone.append(f"{i},{bond[0]},{bond[1]},{a.residue},{b.residue}")
            i += 1

    print("writing backbone bonds")
    with open(f"results/{system}_backbone.csv", "w") as f:
        f.write("\n".join(ls_backbone))

    print("writing peptide bonds")
    with open(f"results/{system}_peptide.csv", "w") as f:
        f.write("\n".join(ls_peptide))

    plumed_path = f"{example}/assets/{system}_plumed.dat"
    plumed_ls = []
    for i, (a, b) in enumerate(backbone):
        plumed_ls.append(f"d{i}: DISTANCE ATOMS={a},{b}")

    plumed_ls.append("")
    plumed_ls.append(
        "PRINT ARG="
        + ",".join([f"d{i}" for i, _ in enumerate(backbone)])
        + f", STRIDE={stride} FILE=distances.dat"
    )

    print("writing plumed.dat")
    with open(plumed_path, "w") as f:
        f.write("\n".join(plumed_ls))


def fill_templates(
    example: Path,
    systems: list[str],
    forces_nN: list[int]|list[float],
    phs: list[float],
    rate_types: list[str],
    shears: list[str] = [""],
    use_cluster: bool = False,
):
    md_template = Template(Path(f"{example}/assets/md.template.mdp").read_text())
    kimmdy_template_eq = Template(
        Path(f"{example}/assets/kimmdy_just_eq.template.yml").read_text()
    )
    kimmdy_template_react = Template(
        Path(f"{example}/assets/kimmdy_just_reactions.template.yml").read_text()
    )
    kimmdy_template_relax = Template(
        Path(f"{example}/assets/kimmdy_just_reactions_relax.template.yml").read_text()
    )

    lincs_section = """
    constraint_algorithm    = lincs
    constraints             = h-bonds
    lincs_iter              = 2
    lincs_order             = 4
    """.strip()
    no_lincs_section = """
    constraints             = none
    """.strip()

    xtc_out = 100
    nst_relax = 10000
    nst_pull = 1000000  # 2 ns
    slow_growth_section = f"""
    free-energy = yes
    init-lambda = 0.0
    delta-lambda = {1/nst_relax}
    """.strip()

    for system in systems:
        for force in forces_nN:
            if system == "triple":
                force_mult = 3
            elif system == "collagen":
                force_mult = 33 * 3  # 33 triplehelices in the fibril
            else:
                force_mult = 1
            force_gmx = 602 * force_mult * force  # force in gromacs units kJ/mol/nm
            morse = "yes"
            if force >= 2 and system == "collagen":
                morse = "no"
            nst_eq = 1000000  # 2 ns
            eq = md_template.safe_substitute(
                {
                    "force": force_gmx,
                    "morse": morse,
                    "dt": 0.002,
                    "nsteps": nst_eq,
                    "xtc_out": xtc_out,
                    "additional": "",
                    "lincs": lincs_section,
                    "temperature": T,
                }
            )
            with open(f"{example}/eq_{system}_{force}nN.mdp", "w") as f:
                f.write(eq)

            pull = md_template.safe_substitute(
                {
                    "force": force_gmx,
                    "morse": morse,
                    "dt": 0.002,
                    "nsteps": nst_pull,
                    "xtc_out": xtc_out,
                    "additional": "",
                    "lincs": lincs_section,
                    "temperature": T,
                }
            )
            with open(f"{example}/pull_{system}_{force}nN.mdp", "w") as f:
                f.write(pull)

            if system == "collagen" and force != 0:
                md_shear_template = Template(
                    Path(f"{example}/assets/md-shear.template.mdp").read_text()
                )
                nst_eq_shear = 1000000  # 2 ns
                if force > 0.6 and force < 1.0 and system == "collagen":
                    nst_eq_shear = 1000000  # 2 ns
                eq_shear = md_shear_template.safe_substitute(
                    {
                        "force": force_gmx,
                        "morse": morse,
                        "dt": 0.002,
                        "nsteps": nst_eq_shear,
                        "xtc_out": xtc_out,
                        "additional": "",
                        "lincs": lincs_section,
                        "temperature": T,
                    }
                )
                pull_shear = md_shear_template.safe_substitute(
                    {
                        "force": force_gmx,
                        "morse": morse,
                        "dt": 0.002,
                        "nsteps": nst_pull,
                        "xtc_out": xtc_out,
                        "additional": "",
                        "lincs": lincs_section,
                        "temperature": T,
                    }
                )
                # similar to
                # /hits/fast/mbm/rennekbt/BC_newEdis_shear/create_shearpulling_v14.py
                average_force = 602 * 3 * force
                std_dev = average_force / 3.0
                random.seed(4234)
                n_groups_per_side = 33

                # generate balanced random forces
                # pull-group1-name         = ACE_0
                # pull-group2-name         = NME_0
                # pull-group3-name         = ACE_1
                # pull-group4-name         = NME_1

                # pull-coord1-k            = -K_PLACEHOLDER
                # pull-coord1-groups       = 0 2
                # pull-coord2-k            = K_PLACEHOLDER
                # pull-coord2-groups       = 0 1
                # pull-coord3-k            = -K_PLACEHOLDER

                ace = []
                nme = []
                for _ in range(n_groups_per_side):
                    ace.append(random.gauss(average_force, std_dev))
                    nme.append(random.gauss(average_force, std_dev))

                # check for negative forces
                for i, f in enumerate(ace):
                    if f < 0:
                        print(
                            f"warning: negative random pulling strength. Replaced with 0"
                        )
                        ace[i] = 0
                for i, f in enumerate(nme):
                    if f < 0:
                        print(
                            f"warning: negative random pulling strength. Replaced with 0"
                        )
                        nme[i] = 0

                # balance with groups from the middle of the pack
                sum_ace = sum(ace)
                sum_nme = sum(nme)
                i = 13
                delta = sum_ace - sum_nme
                if delta / 2 > ace[i] or delta / 2 > nme[i]:
                    print(
                        "warning: Cannot balance pulling sides since delta is higher than last pull group force"
                    )
                ace[i] = ace[i] - delta / 2
                nme[i] = nme[i] + delta / 2

                # replace placeholders in mdp file strings
                for name, mdp_string in [("eq", eq_shear), ("pull", pull_shear)]:
                    i = 0
                    ls = []
                    a = 0
                    n = 0
                    for l in mdp_string.splitlines():
                        if "K_PLACEHOLDER" in l:
                            if i % 2 == 0:
                                assert "-K_P" in l
                                l = l.replace("K_PLACEHOLDER", str(int(round(ace[a]))))
                                a += 1
                            else:
                                assert " K_P" in l
                                l = l.replace("K_PLACEHOLDER", str(int(round(nme[n]))))
                                n += 1
                            i += 1
                        ls.append(l)
                    mdp_string = "\n".join(ls)
                    with open(
                        f"{example}/{name}_{system}_{force}nN_shear.mdp", "w"
                    ) as f:
                        f.write(mdp_string)

            relax = md_template.safe_substitute(
                {
                    "force": force_gmx,
                    "morse": "no",
                    "dt": "0.0005",
                    "nsteps": nst_relax,
                    "xtc_out": xtc_out,
                    "lincs": no_lincs_section,
                    "additional": slow_growth_section,
                    "temperature": T,
                }
            )
            with open(f"{example}/relax_{system}_{force}nN.mdp", "w") as f:
                f.write(relax)

            if system == "collagen":
                if force <= 0.8:
                    # nvt of collagen is at 0 nN
                    gro = "nvt"
                else:
                    # npt is already equilibrated at 1 nN
                    gro = "npt"
            else:
                gro = f"{system}_npt"
            if use_cluster:
                # PLUMED crashes with more than 1 ntmpi
                md_flags = "-maxh 23 -dlb yes -ntmpi 1 -nt 40"
                max_hours = 23
            else:
                md_flags = "-ntmpi 1 -nt 10 -dlb yes --gpu_id 0"
                max_hours = 0

            for shear in shears:
                if not (system == "collagen" and force != 0):
                    continue
                yml_just_eq = kimmdy_template_eq.safe_substitute(
                    {
                        "system": system,
                        "force": force,
                        "gro_start": gro,
                        "md_flags": md_flags,
                        "max_hours": max_hours,
                        "shear": shear,
                    }
                )
                with open(
                    f"{example}/kimmdy_eq_{system}_{force}nN{shear}.yml", "w"
                ) as f:
                    f.write(yml_just_eq)

            for ph in phs:
                for rate_type in rate_types:
                    for shear in shears:
                        yml_just_reactions = kimmdy_template_react.safe_substitute(
                            {
                                "system": system,
                                "force": force,
                                "external_force": -1,
                                "ph": ph,
                                "temperature": T,
                                "max_hours": max_hours,
                                "rate_type": rate_type,
                                "use_theoretical_rates": rate_type == "theo",
                                "bondstats_at_0": f"./assets/{system}_bondstats.csv",
                                "shear": shear,
                            }
                        )
                        with open(
                            f"{example}/kimmdy_just_reactions_{system}_{force}nN_{ph}pH_{rate_type}{shear}.yml",
                            "w",
                        ) as f:
                            f.write(yml_just_reactions)

                        yml_relax = kimmdy_template_relax.safe_substitute(
                            {
                                "system": system,
                                "force": force,
                                "external_force": -1,
                                "ph": ph,
                                "temperature": T,
                                "md_flags": md_flags,
                                "max_hours": max_hours,
                                # TODO: apply to template:
                                "bondstats_at_0": f"./assets/{system}_bondstats.csv",
                                "shear": shear,
                            }
                        )
                        with open(
                            f"{example}/kimmdy_relax_{system}_{force}nN_{ph}pH{shear}.yml",
                            "w",
                        ) as f:
                            f.write(yml_relax)


def fill_templates_collagen(example: Path, forces_nN: list[int], phs: list[float]):
    md_template_path = Path(f"{example}/assets/md.template.mdp")
    md_template = Template(md_template_path.read_text())

    kimmdy_template_just_eq = Template(
        Path(f"{example}/assets/kimmdy-just-eq.template.yml").read_text()
    )
    kimmdy_template_just_reactions = Template(
        Path(f"{example}/assets/kimmdy-just-reactions.template.yml").read_text()
    )
    kimmdy_template_relax = Template(
        Path(f"{example}/assets/kimmdy-just-reactions-relax.template.yml").read_text()
    )

    lincs_section = """
    constraint_algorithm    = lincs
    constraints             = h-bonds
    lincs_iter              = 2
    lincs_order             = 4
    """.strip()
    no_lincs_section = """
    constraints             = none 
    """.strip()

    xtc_out = STRIDE
    nst_relax = 10000
    nst_pull = 1000000
    slow_growth_section = f"""
    free-energy = yes
    init-lambda = 0.0
    delta-lambda = {1/nst_relax}
    """.strip()

    for force in [0, 1]:
        force_gmx = force * 59598
        eq = md_template.substitute(
            {
                "force": force_gmx,
                "morse": "yes",
                "dt": 0.002,
                "nsteps": 1000000,  # 2 ns
                "xtc_out": xtc_out,
                "additional": "",
                "lincs": lincs_section,
            }
        )
        with open(f"{example}/eq-{force}nN.mdp", "w") as f:
            f.write(eq)

        pull = md_template.substitute(
            {
                "force": force_gmx,
                "morse": "yes",
                "dt": "0.002",
                "nsteps": nst_pull,
                "xtc_out": xtc_out,
                "additional": "",
                "lincs": lincs_section,
            }
        )
        pull_path = f"{example}/pull-{force}nN.mdp"
        with open(pull_path, "w") as f:
            f.write(pull)

        relax = md_template.safe_substitute(
            {
                "force": force_gmx,
                "morse": "no",
                "dt": "0.0005",
                "nsteps": nst_relax,
                "xtc_out": xtc_out,
                "lincs": no_lincs_section,
                "additional": slow_growth_section,
            }
        )
        relax_path = f"{example}/relax-{force}nN.mdp"
        with open(relax_path, "w") as f:
            f.write(relax)

        if force == 0:
            # nvt is at 0 nN
            gro = "nvt"
        else:
            # npt is already equilibrated at 1 nN
            gro = "npt"

        yml_just_eq = kimmdy_template_just_eq.safe_substitute(
            {"name": f"run_eq-{force}nN", "force": force, "gro_start": gro}
        )
        yml_path_just_eq = f"{example}/kimmdy_eq-{force}nN.yml"
        with open(yml_path_just_eq, "w") as f:
            f.write(yml_just_eq)

        yml_just_reactions = kimmdy_template_just_reactions.safe_substitute(
            {
                "name": f"run_{force}nN",
                "force": force,
                "external_force": -1,
                "ph": 7.4,
                "eq_bond_lengths": f"./assets/collagen_bond_lengths.csv",
            }
        )
        yml_path = f"{example}/kimmdy_reactions-{force}nN.yml"
        with open(yml_path, "w") as f:
            f.write(yml_just_reactions)


def gather_0_dists_and_forces(example: Path, system: str, dt: float = 0):
    run = f"{example}/run_eq_{system}_0nN"
    top_name = f"{system}.top"
    print(f"Reading bondstats at 0 nN of {system}")
    plumed_in = Path(f"{run}/0_setup/{system}_plumed.dat")
    plumed_out = Path(f"{run}/2_pull/distances.dat")
    edissoc_dat = Path(f"{example}/assets/edissoc.dat")
    top = Topology.from_path(f"{run}/0_setup/{top_name}")

    stats = calculate_bondstats(
        top=top,
        plumed_in=plumed_in,
        plumed_out=plumed_out,
        dt=dt,
        edissoc_dat=edissoc_dat,
    )
    bondstats_to_csv(stats=stats, path=f"{example}/assets/{system}_bondstats.csv")


def slurm_dispatch(
    example: str | Path,
    job: str,
):
    """
    Dispatch a job to the slurm cluster
    """
    sh(f"ssh {CLUSTER} 'cd {example} && sbatch jobscript-{job}.sh'")
