# example

![Figure from Hamilton Population Genetics book](dna_polymorphism.jpg)

# definitions

## $\theta_w$ Watterson's theta

the number of segregating sites divided by the harmonic number of (n-1) per nucleotide site

* $S$: The number of segregating sites
* $L$: The total number of sites
* $p_S$: The number of segregating sites per nucleotide site
* $n$: The number of lineages sampled (i.e. sample size)

$$
p_S = S / L
$$

$$
\theta_w = \frac{p_S}{\sum_{k=1}^{n-1} \frac{1}{k}}
$$


## $\theta_\pi$ nucleotide diversity

the average pairwise differences in a sample of DNA sequences per nucleotide site

* $n$: The number of lineages sampled (i.e. sample size)
* $d_{ij}$: The number of differences in DNA sequence samples $i$ and $j$
* $C_n^2$: The number of combinations of $n$ samples taken 2 times (number of pairwise DNA sequence comparisons)
* $L$: The total number of sites


$$
\theta_\pi = \frac{\sum^{n-1}_{i=1}\sum^n_{j>i}d_{ij}}{C_n^2L}
$$

To rewrite the above definitions by allele frequencies for **biallelic SNPs**.

* $p_l$: allele $A$'s frequency at SNP position $l$
* $q_l$: allele $a$'s frequency at SNP position $l$
* $np_l$: the count of samples with allele $A$ at SNP position $l$
* $nq_l$: the count of samples with allele $a$ at SNP position $l$
* $p_l + q_l = 1$

$$
\begin{align}
\theta_\pi &= \frac{\sum^{n-1}_{i=1}\sum^n_{j>i}d_{ij}}{\frac{n(n-1)}{2}L}\\
& = \frac{\sum^L_{l=1}(np_l)(nq_l)}{\frac{n(n-1)}{2}L} \\
& = \frac{n}{(n-1)L}*\sum^L_{l=1}2p_lq_l
\end{align}
$$

> In the step(2) here, we changed the thinking from counting pairwise difference for a sequence length of $L$ to counting the pairwise difference for a given SNP location $l$. Now the problem becomes the hypergeometric distribution "draw two balls from a mixture of A and a balls ($n$ balls, number of A balls is $c_A=np_l$), and only one ball is A", the probability is:
> $\frac{C^1_{c_A}C^1_{c_a}}{C^2_n}$

For non-biallelic alleles, when $n$ is large enough ($n \sim n-1$) we can write:

* $M_l$: The number of alleles at SNP position $l$
* $p_{ml}$: allele $A_m$'s frequency at SNP position $l$
* $\sum^{M_l}_{m=1} p_{ml} = p_{1l} + p_{2l} + ... + p_{M_ll} = 1$

$$
\theta_\pi = \frac{1}{L}*\sum^L_{l=1}(1-\sum^{M_l}_{m=1}p^2_{ml})
$$

**Q: is the assumption correct?**

## implications for `extinction-sim`

* `rasterN`: should not be used. Each extinction event is sampling differing number of lineages.
* `L`: is probably underestimated as the number of raster layers are the total number of segregating sites.

# predictions about $\theta_w$ and $\theta_\pi$

Assuming Hardy-Winberg Equilibrium: **Q: maybe a less stringent assumption suffice?**

1. $\theta_\pi = $ the number of heterozygote genotypes (tally `0|1` or `0/1` genotypes) per nucleotide site.



