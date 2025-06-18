from pathlib import Path

import pandas as pd
from kimmdy.recipe import Break, DeferredRecipeSteps, RecipeCollection
from kimmdy_hydrolysis.reaction import read_sasa_free
import polars as pl
import numpy as np


def prepare_data():
    # Using polars instead of pandas due to the size of the data
    rates = pl.scan_csv("./results/rates.csv")

    exp_rates = pd.read_csv("./assets/experimental-rates.csv")
    exp_rates = exp_rates.query("force in [0, 1, 2, 3]")
    exp_rates["rate_type"] = "exp"
    exp_rates["ph"] = 7.4
    exp_rates["system"] = "single peptide"
    exp_rates["reaction"] = "hydrolysis"
    exp_rates["shear"] = False
    exp_rates = exp_rates.rename(columns={"rate": "rate_ps"})
    exp_rates["rate"] = exp_rates["rate_s"]

    highest_rates = (
        rates.filter(pl.col("rate") > 1e-30)
        .group_by(
            [
                "reaction",
                "system",
                "force",
                "ph",
                "rate_type",
                "shear",
            ]
        )
        .agg(pl.max("rate").alias("rate"))
        .collect()
        .to_pandas()
    )
    highest_rates["rounded_rate"] = highest_rates["rate"].apply(
        lambda x: "{:.1e}".format(x)
    )

    np.random.seed(1)

    sampled_rates = (
        rates.filter(pl.col("rate") > 1e-30)
        .group_by(
            [
                "force",
                "ph",
                "system",
                "reaction",
                "rate_type",
                "shear",
            ]
        )
        .map_groups(lambda df: df.sample(min(len(df), 50000)), schema=None)
        .collect()
        .to_pandas()
    )

    for df in [sampled_rates, exp_rates, highest_rates]:
        df["system"] = df["system"].str.title()
        df["reaction_type"] = df["reaction"] + " " + df["rate_type"]
        df["cat"] = df["reaction_type"] + df["shear"].astype(str)
        df["shear"] = df.apply(
            lambda x: True if x["force"] == 0.0 else x["shear"], axis=1
        )

    highest_rates.to_csv(
        "./results/highest_rates.csv",
        index=False,
    )

    exp_rates.to_csv(
        "./results/exp_rates.csv",
        index=False,
    )

    sampled_rates.to_csv(
        "./results/sample_rates.csv",
        index=False,
    )


def get_rates(
    root: str | Path,
    systems: list[str],
    forces: list[int],
    phs: list[float],
    rate_types: list[str],
    shears: list[str] = [""],
) -> None:
    single_bonds = pd.read_csv("./results/single_backbone.csv")
    triple_bonds = pd.read_csv("./results/triple_backbone.csv")
    collagen_bonds = pd.read_csv("./results/collagen_backbone.csv")

    single_peptide_bonds = pd.read_csv("./results/single_peptide.csv")
    triple_peptide_bonds = pd.read_csv("./results/triple_peptide.csv")
    collagen_peptide_bonds = pd.read_csv("./results/collagen_peptide.csv")

    Path("./results/rates").mkdir(parents=True, exist_ok=True)

    def read_rates(path, force, system, ph, rate_type, shear, recompute: bool = False):
        print(f"Reading {path}")

        name = f"{system}_{force}nN_{ph}pH_{rate_type}{shear}"
        path_out = Path(f"./results/rates/{name}.csv")
        if path_out.exists() and not recompute:
            print(f"Already exists {path_out}")
            return

        collection, _ = RecipeCollection.from_csv(path)

        # assuming order of rates of homolysis follows order of backbone bonds
        # and order of hydrolysis follows order of peptide bonds
        i_homolysis = 0
        i_hydrolysis = 0

        if system == "single":
            bonds = single_bonds
            peptide_bonds = single_peptide_bonds
        elif system == "triple":
            bonds = triple_bonds
            peptide_bonds = triple_peptide_bonds
        elif system == "collagen":
            bonds = collagen_bonds
            peptide_bonds = collagen_peptide_bonds
        else:
            raise ValueError(f"Unknown system: {system}")

        rates = []
        for recipe in collection.recipes:
            steps = recipe.recipe_steps
            times = recipe.timespans
            if isinstance(steps, list):
                # homolysis just has Break and Relax as steps
                # and only one timespan and rate per recipe
                # (=per bond)
                # and homolysis does not differ between rate_types
                # (which is just a setting for hydrolysis)
                # so we can skip one or the other rate_types
                if rate_type == "exp":
                    continue
                # homolysis will be rate_type "theo"

                assert len(recipe.rates) == 1
                step = steps[0]
                assert isinstance(step, Break)
                ai = step.atom_id_1
                aj = step.atom_id_2
                row = bonds.iloc[i_homolysis]
                c = row["c"]
                n = row["n"]
                i_homolysis += 1
                time = times[0][1]
                rate = recipe.rates[0]
                rates.append(
                    {
                        "force": force,
                        "system": system,
                        "ph": ph,
                        "reaction": "homolysis",
                        "rate_type": rate_type,
                        "shear": bool(shear),
                        "time": time,
                        "rate": rate,
                        "ai": ai,
                        "aj": aj,
                        "c": c,
                        "n": n,
                    }
                )

            if isinstance(steps, DeferredRecipeSteps):
                # hydrolysis uses the DeferredRecipeSteps
                # and has multiple rates and times per recipe
                # (=per peptide bond)
                steps.key = int(steps.key)
                row = peptide_bonds.iloc[i_hydrolysis]
                i_hydrolysis += 1
                ai = row["ai"]
                aj = row["aj"]
                c = row["c"]
                n = row["n"]
                for rate, time in zip(recipe.rates, recipe.timespans):
                    rates.append(
                        {
                            "force": force,
                            "system": system,
                            "ph": ph,
                            "reaction": "hydrolysis",
                            "rate_type": rate_type,
                            "shear": bool(shear),
                            "time": time[1],
                            "rate": rate,
                            "ai": ai,
                            "aj": aj,
                            "c": c,
                            "n": n,
                        }
                    )

        # write inidividual rates to csvs
        rates = pd.DataFrame(rates)
        rates.to_csv(path_out, index=False)

    for system in systems:
        for force in forces:
            for ph in phs:
                for rate_type in rate_types:
                    for shear in shears:
                        if system == "collagen":
                            example = f"{root}/examples/collagen-hydrolysis"
                        else:
                            example = f"{root}/examples/triplehelix-hydrolysis"
                        path = Path(
                            f"{example}/run_react_{system}_{force}nN_{ph}pH_{rate_type}{shear}/3_decide_recipe/recipes.csv"
                        )
                        if not path.exists():
                            continue
                        read_rates(
                            path=path,
                            force=force,
                            system=system,
                            ph=ph,
                            rate_type=rate_type,
                            shear=shear,
                        )

    # combine all rates into one dataframe
    rates = pd.DataFrame()
    for system in systems:
        for force in forces:
            for ph in phs:
                for rate_type in rate_types:
                    for shear in shears:
                        name = f"{system}_{force}nN_{ph}pH_{rate_type}{shear}"
                        path = Path(f"./results/rates/{name}.csv")
                        if not path.exists():
                            continue
                        rates = pd.concat([rates, pd.read_csv(path)], ignore_index=True)

    rates["rate"] *= 2e12  # convert from 1/ps to 1/s
    rates["reaction_type"] = (
        rates["reaction"] + " " + rates["rate_type"] + " " + rates["shear"].astype(str)
    )
    rates["system"] = rates["system"].replace("collagen", "fibril")
    rates["system"] = (
        rates["system"]
        .replace("single", "single peptide")
        .replace("triple", "triple helix")
    )
    rates.to_csv("./results/rates.csv", index=False)
    # don't return the rates, such that it can be garbage collected


def get_sasa(
    root: str | Path,
    systems: list[str],
    forces: list[str | int | float],
    shears: list[str] = [""],
):
    ls = []

    def read_sasa(path, force, system, shear):
        result = read_sasa_free(path)
        if result is None:
            print(f"Error reading {path}")
            return
        meta, times, sasas = result

        for i, time in enumerate(times):
            for k in sasas.keys():
                sasa = sasas[k]
                ls.append(
                    {
                        "system": system,
                        "force": force,
                        "shear": bool(shear),
                        "time": time,
                        "ai": k,
                        "sasa": sasas[k][i],
                    }
                )

    for system in systems:
        for force in forces:
            for shear in shears:
                if system == "collagen":
                    example = f"{root}/examples/collagen-hydrolysis"
                else:
                    example = f"{root}/examples/triplehelix-hydrolysis"
                path = Path(
                    f"{example}/run_eq_{system}_{force}nN{shear}/2_pull/.kimmdy.sasa"
                )
                if not path.exists():
                    continue
                print(f"Reading {path}")
                read_sasa(path=path, force=force, system=system, shear=shear)

    df = pd.DataFrame(ls)
    df.to_csv("./results/sasa.csv", index=False)
    return df
