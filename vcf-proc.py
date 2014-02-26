#!/software/bin/python
# Aylwyn Scally 2013
#TODO: handle ./. output e.g. from vcf-merge

import sys
import optparse
#import os
import re
#import os.path
from numpy import median
import logging
from logging import error, warning, info, debug, critical

loglevel = logging.WARNING

# defaults
VCF_FIXEDCOLS = 9
#HWEMIN = 0.05

#	fout.write('\t'.join([str(x) for x in args]) + '\n')

def genotype(al1, al2):
	return(''.join(sorted([al1, al2])))

def hapsitetypes(n):
# haploid allele combinations
	for gt in ['0', '1']:
		if n > 1:
			for cc in hapsitetypes(n-1):
				yield(''.join((gt, cc)))
		else:
			yield(gt)


def dipsitetypes(n):
# diploid genotype combinations
# to get segregating sites drop first and last elements of returned vals
	for gt in ['00', '01', '11']:
		if n > 1:
			for cc in dipsitetypes(n-1):
				yield('-'.join((gt, cc)))
		else:
			yield(gt)

def phases(gt):
# diploid allele combinations
	yield(gt)
	if gt[0] != gt[1]:
		yield(gt[::-1])

def phasecombs(gts, sep = ''):
# different phase combinations of diploid genotypes in gts 
	if len(gts) > 1:
		for cc in phasecombs(gts[1:], sep):
			for ph in phases(gts[0]):
				yield(sep.join((ph + cc)))
	else:
		for ph in phases(gts[0]):
			yield(ph)

#test
tgt = ['AB', 'CD', 'EF']
assert ','.join(phasecombs(tgt)) == 'ABCDEF,BACDEF,ABDCEF,BADCEF,ABCDFE,BACDFE,ABDCFE,BADCFE'

#p = optparse.OptionParser(usage = '%prog [file] [--segsites] [--diversity] [-c cond_sample]')
p = optparse.OptionParser()
p.add_option('-H', '--header', action='store_true', default = False)
p.add_option('--depths', action='store_true', default = False, help = 'output per-individual read depths')
p.add_option('--indels', action='store_true', default = False, help = 'include indels')
p.add_option('--segsites', action='store_true', default = False, help = 'output segregating site flag')
p.add_option('--mediandep', action='store_true', default = False, help = 'output median read depth')
p.add_option('--diversity', action='store_true', default = False, help = 'output diversity statistics')
p.add_option('-c', '--condsamp', default = '', help = 'condition on het in sample CONDSAMP')
p.add_option('--rates', action='store_true', default = False, help = 'output per-individual het and hom-non-ref rates')
p.add_option('--ratesonly', action='store_true', default = False, help = 'as --rates; suppress per-site output')
p.add_option('--vars', action='store_true', default = False, help = 'output variant sites only')
p.add_option('--segsep', action='store_true', default = False, help = 'output seg site separations (e.g. for input to msmc)')
p.add_option('--replacecalls', default='', help = 'file of replacement records (e.g. phased SNP calls). NOTE: no checking is done to ensure samples match')
p.add_option('--alleles', action='store_true', default = False, help = 'output alleles')

opt, args = p.parse_args()
if opt.ratesonly:
	opt.rates = True

if args:
	fin = open(args[0])
else:
	fin = sys.stdin
fout = sys.stdout

if opt.replacecalls:
	reprecord = {}
	for line in open(opt.replacecalls):
		if line.startswith('#'):
			continue
		tok = line.split()
		reprecord[tok[0]+tok[1]] = tok[VCF_FIXEDCOLS:]


for line in fin:
#	if not opt.noheader:
#		sys.stdout.write(line)
	if line.startswith('#CHROM'):
		tok = line.split()
		nsamp = len(tok) - VCF_FIXEDCOLS
		sampnames = tok[VCF_FIXEDCOLS:]
		if opt.condsamp:
			for i, name in enumerate(sampnames):
				if name.find(condsamp) >= 0:
					condcol = i + VCF_FIXEDCOLS
#				   print(i, condsamp, tok[condcol])
					break
		break

if opt.header:
	fout.write('\t'.join(['#SAMPLES'] + sampnames) + '\n')

if opt.rates:
	hetcount = [0 for i in range(nsamp)]
	hnrcount = [0 for i in range(nsamp)]
	depcount = [0 for i in range(nsamp)]
	nsites = 0
	firstcoord = ''

lastchr = ''
for line in fin:
	if line.startswith('#'):
#		sys.stdout.write(line)
		continue
	#TODO: handle new chr: reset sep, lastpos
	if tok[0] != lastchr:
		sep = 0
		lastpos = -1
		lastchr = tok[0]

	tok = line.split()

	if opt.replacecalls and tok[0]+tok[1] in reprecord:
		tok[VCF_FIXEDCOLS:] = reprecord[tok[0]+tok[1]]
	
	if opt.vars and tok[4] == '.':
		continue

#	sys.stdout.write(line)
	if not opt.indels and (len(tok[3]) > 1 or re.search('INDEL', tok[7])):
		continue

#	if len(tok[4]) > 1: # multiallelic sites
#		continue

#	if nsamp > 1:
#		if opt.HWEmin > 0.0:
#			hwep = re.search('HWE=([0-9.e\-]+)', tok[7])

	if opt.condsamp and not tok[condcol].startswith('0/1'): # condition on variant in cond_sample
		continue

	if opt.depths or opt.depths or opt.mediandep:
		formats = tok[VCF_FIXEDCOLS - 1].split(':')
		depfield = formats.index('DP')
		dep = [x.split(':')[depfield] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]

	if opt.depths:
		fout.write('\t'.join(tok[0:2] + dep) + '\n')
		continue

	outvals = tok[0:2]

	if (tok[4] == "."):
		varvals = ['00' for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
		segsite = False
	else:
		varvals = [x[0] + x[2] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
##		varvals = [x.split(':')[0].replace('/', '') for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
#		if sdafilt:
#			if (varvals[-1] == '01' and varvals.count('00') < nfsamp) or (varvals[-1] == '11' and varvals.count('01')): #exclude lines with sda in last sample
#				continue
#	outvals.append('-'.join(varvals[:nfsamp]))
		segsite = '01' in varvals

	if opt.alleles: #TODO
		sitealleles = [tok[3]] + tok[4].split(',')
		varalleles = [sitealleles[int(x[0])] + sitealleles[int(x[1])] for x in varvals]  
		varvals = varalleles  
#		print([sitealleles, varvals])

	if opt.segsep:
		pos = int(tok[1])
		if pos > lastpos:
			sep += 1 # assumes all callable sites are in vcf TODO otherwise sep = pos - lastpos
		lastpos = pos
#	ideally use FQ > 0 but vcf-subset doesn't currently adjust this
#		fq = re.search('FQ=(-?[0-9.]+)', tok[7])
#		if fq > 0: #het
#		sys.stdout.write(line)
		if segsite:
			if len(varvals) > 1: #
				allstr = ','.join(phasecombs(varvals))
				fout.write('\t'.join(outvals + [str(sep)] + [allstr]) + '\n')
			else:
				fout.write('\t'.join(outvals + [str(sep)] + varvals) + '\n')
#			outvals.append(str(sep))
			sep = 0
		continue

	outvals.append('-'.join(varvals))

	if opt.rates:
		hetcount = [hetcount[i] + (varvals[i] == '01') for i in range(nsamp)]
		hnrcount = [hnrcount[i] + (varvals[i] == '11') for i in range(nsamp)]
		depcount = [depcount[i] + int(dep[i]) for i in range(nsamp)]
		nsites += 1
		if not firstcoord:
			firstcoord = ':'.join(tok[0:2])
		lastcoord = ':'.join(tok[0:2])

	if not opt.ratesonly:
		if opt.mediandep:
			outvals.append('%d' % median(dep))

		if opt.segsites:
#			outvals.append(str(depfield - 1))
			outvals.append(str(int(segsite)))

		if opt.diversity:
			alfreq = eval(re.search('AF1=([0-9.e\-]+)', tok[7]).group(1))# eval since vcf uses sci notation
			alcount = int(re.search('AC1=([0-9]+)', tok[7]).group(1))
			neidiv = 1 - alfreq**2 - (1 - alfreq)**2
			outvals += [str(nsamp), str(alfreq), str(neidiv), str(alcount), str(alcount/float(nsamp))]

		fout.write('\t'.join(outvals) + '\n')

if opt.rates:
	fout.write('\t'.join(['#REGION'] + [firstcoord + '-' + lastcoord]) + '\n')
	fout.write('\t'.join(['#NSITES'] + [str(nsites)]) + '\n')
	fout.write('\t'.join(['#MEANDEP'] + [str(float(x)/nsites) for x in depcount]) + '\n')
	fout.write('\t'.join(['#HETCOUNT'] + [str(x) for x in hetcount]) + '\n')
	fout.write('\t'.join(['#HETRATE'] + [str(float(x)/nsites) for x in hetcount]) + '\n')
	fout.write('\t'.join(['#HNRCOUNT'] + [str(x) for x in hnrcount]) + '\n')
	fout.write('\t'.join(['#HNRRATE'] + [str(float(x)/nsites) for x in hnrcount]) + '\n')
	fout.write('\t'.join(['#HNR/HET'] + [str(float(hnrcount[i])/hetcount[i]) for i in range(nsamp)]) + '\n')
#	fout.write('\t'.join(outvals) + '\n')
