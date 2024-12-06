#!/usr/bin/env python3

# Based on: https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/h5m2csv.py
# h5m2csv.py: Convert HDF5 SNP matrix to CSV
#
# (c) 2021 by Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de),
# Max Planck Institute for Developmental Biology,
# Tuebingen, Germany
# Updated by meixilin to format inputs for mar package
# 2024-11-29

import h5py, numpy, pandas

f = h5py.File('1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')

# Check sample id and cleanup accessions
lonlat = pandas.read_csv('1001_accessions.csv', delimiter=',', usecols=[0,3,5,6], quotechar='"', header=None)
lonlat.columns = ['ID','CNTY','LAT', 'LON']
# Check that the csv file is the same as the hd5 file
accessions = f['accessions'][:].astype(int)
if not numpy.array_equal(lonlat['ID'].values, accessions):
	raise ValueError("Accession numbers in csv and hdf5 file do not match")
# Remove samples with missing coordinates and samples in USA, CAN and JPN
lonlatbad = lonlat[lonlat.isna().any(axis=1) | lonlat['CNTY'].isin(['USA', 'CAN', 'JPN'])]
print(lonlatbad) # 131 samples
# badsampidx = (array([439, 499, 501, 503]),)
badsampidx = numpy.where(numpy.isin(accessions, lonlatbad['ID'].values))
accessions = numpy.delete(accessions, badsampidx)
# Write sample ids
numpy.savetxt('1001g_accessions.txt', accessions, fmt='%d')
# Write lonlat data
lonlatgood = lonlat[~lonlat['ID'].isin(lonlatbad['ID'])][['ID', 'LON', 'LAT']]
if not numpy.array_equal(lonlatgood['ID'].values, accessions):
	raise ValueError("Accession numbers in csv and hdf5 file do not match")
lonlatgood.to_csv('1001g_lonlat.txt', index=False, header=True, sep='\t')

# Get all SNP positions for all chromosomes (len=10709949)
positions = f['positions'][:]
# Array of tupels with start/stop indices for each chromosome
chr_regions = f['positions'].attrs['chr_regions']
chr1end = chr_regions[0][1] # 2597825

# Only sample 10000 SNPs from the first chromosome
nsnps = 10000
randidx = numpy.random.choice(chr1end, size = nsnps, replace = False)
# Check the SNPs that are homozygous after removing the samples with missing coordinates
for ii, idx in enumerate(randidx):
	snps = numpy.delete(f['snps'][idx], badsampidx)
	if numpy.sum(snps) == 0 or numpy.sum(snps) == len(snps):
		while numpy.sum(snps) == 0 or numpy.sum(snps) == len(snps):
			print((ii, idx))
			notrandidx = numpy.setdiff1d(numpy.arange(chr1end), randidx)
			idx = numpy.random.choice(notrandidx, size=1, replace=False)[0]
			snps = numpy.delete(f['snps'][idx], badsampidx)
		randidx[ii] = idx
randidx = numpy.sort(randidx)
if len(numpy.unique(randidx)) != nsnps:
	raise ValueError("Number of SNPs is not correct")
numpy.savetxt('1001g_chr1.10k_randidx.txt', randidx, fmt='%d')

# Write chromosome positions in the sampled 10000 SNPs in chr1
with (open('1001g_chrpos.txt', 'w') as fpos):
	fpos.write("CHR\tPOS\n")
	for idx in randidx:
		out = "1\t{0}\n".format(f['positions'][idx])
		fpos.write(out)

# Check genotypes
randsnps = f['snps'][randidx]
print(numpy.unique(randsnps, return_counts = True))

# Remove the SNPs belonging to the samples with missing coordinates
randsnps = numpy.delete(randsnps, badsampidx, axis=1)
# Check the SNPs that are homozygous after removing the samples with missing coordinates
randac = numpy.sum(randsnps, axis=1)
print(numpy.unique(randac, return_counts = True))
if numpy.any(randac == 0) or numpy.any(randac == randsnps.shape[1]):
	raise ValueError("Some SNPs are homozygous")
# Convert 1 to 2 as A. thaliana is diploid homozygous
randsnps[randsnps == 1] = 2
print(randsnps.shape)
numpy.savetxt('1001g_genotypes.txt', randsnps, fmt='%d', delimiter='\t')
