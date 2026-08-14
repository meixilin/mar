#!/usr/bin/env python3
"""Simulate a diploid population under exponential growth with msprime and write MAR-ready inputs.

Writes two files usable as `mar::vcf_parser()` / `mar::lonlat_parser()` inputs:
  - gmexp.vcf.gz     : gzipped VCF of simulated biallelic SNPs
  - gmexp_lonlat.csv : random lon/lat coordinates, header "ID,LON,LAT"

Based on: ARCHIVE/mar_theory/scripts/neutral.ipynb
"""

import gzip

import msprime
import numpy as np

# simulation parameters --------
sample_size = 50  # number of diploid individuals to sample
sequence_length = 5e6
recombination_rate = 1e-8
mutation_rate = 1e-8
initial_size = 10000  # population size at present (time 0)
growth_rate = 0.01  # exponential growth rate per generation (see msprime docs)
end_time = 500  # generations ago when growth stops and ancestral size becomes constant
seed = 1
out_prefix = "gmexp"

# simulate ancestry and mutations --------
# msprime growth_rate: N(t) = initial_size * exp(-growth_rate * t), t = generations before present.
# growth_rate > 0 means the population grew exponentially towards the present.
demography = msprime.Demography()
demography.add_population(name="pop", initial_size=initial_size, growth_rate=growth_rate)
demography.add_population_parameters_change(time=end_time, population="pop", growth_rate=0)

ts = msprime.sim_ancestry(
    samples=sample_size,
    demography=demography,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    random_seed=seed,
)
mts = msprime.sim_mutations(ts, rate=mutation_rate, model="binary", random_seed=seed)

# write vcf --------
vcf_path = f"{out_prefix}.vcf.gz"
with gzip.open(vcf_path, "wt") as f:
    mts.write_vcf(f)

# write lonlat --------
lonlat_path = f"{out_prefix}_lonlat.csv"
rng = np.random.default_rng(seed)
coords = rng.uniform(0, 1000, size=(mts.num_individuals, 2))
with open(lonlat_path, "w") as f:
    f.write("ID,LON,LAT\n")
    for i, (lon, lat) in enumerate(coords):
        f.write(f"tsk_{i},{lon},{lat}\n")

print(f"num_trees: {mts.num_trees}")
print(f"segregating_sites: {mts.segregating_sites(span_normalise = False)}")
print(f"wrote {vcf_path} ({mts.num_individuals} individuals, {mts.num_sites} sites)")
print(f"wrote {lonlat_path}")
