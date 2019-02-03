#!/usr/bin/env python
# Aylwyn Scally 2012-2018
#TODO: handle ./. output e.g. from vcf-merge

import sys
import optparse
#import os
import re
import os.path
from numpy import median
import logging
from logging import error, warning, info, debug, critical
import gzip
from string import maketrans
#import difflib

loglevel = logging.INFO

# defaults
VCF_FIXEDCOLS = 9
PSMC_BINSIZE = 100
PSMC_N_RATIO = 0.9

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
	yield(gt[0] + gt[2])
	if gt[1] != '|' and gt[0] != gt[2]:
#		yield(gt[::-1])
		yield(gt[2] + gt[0])

def phasecombs(gts, sep = ''):
# different phase combinations of diploid genotypes in gts 
	if len(gts) > 1:
		for cc in phasecombs(gts[1:], sep):
			for ph in phases(gts[0]):
				yield(sep.join((ph + cc)))
	else:
		for ph in phases(gts[0]):
			yield(ph)

class FastaStream(object):
	def __init__(self, fh):
		self.fh = fh
		self.ccount = 0

	def write(self, c):
		self.fh.write(c)
		self.ccount += len(c)
		if self.ccount >= 60:
			self.fh.write('\n')
			self.ccount = 0

	def newrec(self, name):
		self.fh.write('\n>%s\n' % str(name))
		self.ccount = 0


#test
tgt = ['A/B', 'C|D', 'E/F']
assert ','.join(phasecombs(tgt)) == 'ABCDEF,BACDEF,ABCDFE,BACDFE'
tgt = ['A/B', 'C/D', 'E/F']
assert ','.join(phasecombs(tgt)) == 'ABCDEF,BACDEF,ABDCEF,BADCEF,ABCDFE,BACDFE,ABDCFE,BADCFE'

#p = optparse.OptionParser(usage = '%prog [file] [--segsites] [--diversity] [-c cond_sample]')
p = optparse.OptionParser()
p.add_option('-H', '--header', action='store_true', default = False)
p.add_option('--debug', action='store_true', default = False)
p.add_option('--depths', action='store_true', default = False, help = 'include read depths in rates')
p.add_option('--depthsonly', action='store_true', default = False, help = 'output per-individual read depths only')
p.add_option('--no_calls', action='store_true', default = False, help = 'suppress output of calls')
p.add_option('--indels', action='store_true', default = False, help = 'include indels')
p.add_option('--segsites', action='store_true', default = False, help = 'output segregating site flag')
p.add_option('--mediandep', action='store_true', default = False, help = 'output median read depth')
p.add_option('--diversity', action='store_true', default = False, help = 'output diversity statistics')
p.add_option('--Fst', action='store_true', default = False, help = 'output Fst (Weir & Cockerham)')
p.add_option('-c', '--condsamp', default = '', help = 'condition on het in sample CONDSAMP')
p.add_option('--rates', action='store_true', default = False, help = 'output per-individual het and hom-non-ref rates')
p.add_option('--ratesonly', action='store_true', default = False, help = 'as --rates; suppress per-site output')
p.add_option('--vars', action='store_true', default = False, help = 'output variant sites only')
p.add_option('--segsep', action='store_true', default = False, help = 'output seg site separations (e.g. for input to msmc)')
p.add_option('--psmcfa', action='store_true', default = False, help = 'output psmcfa (for input to psmc; NOTE psmc only works for two chrs)')
p.add_option('--replacecalls', default='', help = 'vcf.gz file of replacement records (e.g. phased SNP calls). NOTE: no checking is done to ensure samples match')
p.add_option('--callmask', default='', help = 'bed.gz file of uncallable regions')
p.add_option('--pops', default='', help = 'file of sample population assignments (line format: sample pop)')
p.add_option('--aims', action='store_true', default = False, help = 'output AIM sites (requires --pops)')
p.add_option('--alleles', action='store_true', default = False, help = 'output alleles')
p.add_option('--pseudodip', action='store_true', default = False, help = 'create pseudodiploids from consecutive pairs of input samples (assumed haploid so exclude hets)')
p.add_option('--haploid', action='store_true', default = False, help = 'haploid input calls')

opt, args = p.parse_args()

if opt.depthsonly:
	opt.depths = True

if opt.ratesonly:
	opt.rates = True

if opt.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if args:
	fin = open(args[0])
else:
	fin = sys.stdin
fout = sys.stdout

if opt.replacecalls:
	info('Reading replacement calls from %s' % opt.replacecalls)
	reprecord = {}
	for line in gzip.open(opt.replacecalls):
		if line.startswith('#'):
			continue
		tok = line.split()
		reprecord[tok[0]+tok[1]] = tok[VCF_FIXEDCOLS-1:]

if opt.callmask:
	info('Reading uncallable regions from %s' % opt.callmask)
	callmask = {}
	if opt.callmask.endswith('.gz'):
		cmf = gzip.open(opt.callmask)
	else:
		cmf = open(opt.callmask)
	for line in cmf:
		if line.startswith('#'):
			continue
		tok = line.split()
		callmask[int(tok[1]) + 1] = int(tok[2])

if opt.pops:
	info('Reading sample population assignments from %s' % opt.pops)
	samppop = {}
	for line in open(opt.pops):
		if line.startswith('#'):
			continue
		tok = line.split()
		samppop[tok[0]] = tok[1]

	if opt.aims:
		opt.vars = True
		fout.write('\t'.join(['#chr', 'pos'] + poplist + ['AIM_pop']) + '\n')

#	for line in fin:
#		if line.startswith('#'):
#			sys.stdout.write(line)
#			continue
#		tok = line.split()
##		print(tok[0]+tok[1])
#		if tok[0]+tok[1] in reprecord:
##			print(tok[0]+tok[1])
#			tok[VCF_FIXEDCOLS-1:] = reprecord[tok[0]+tok[1]]
#		sys.stdout.write('\t'.join(tok) + '\n')
#	sys.exit(0)

for line in fin:
#	if not opt.noheader:
#		sys.stdout.write(line)
	if line.startswith('#CHROM'):
		tok = line.split()
		nsamp = len(tok) - VCF_FIXEDCOLS
		sampnames = tok[VCF_FIXEDCOLS:]
		sampn = dict([(x[1], int(x[0])) for x in enumerate(sampnames)])

		popsampnums = {}
		if opt.pops: # fill poplist with pops actually present
			poplist = []
			for samp in sampnames:
				if not samp in samppop: # look for matching prefix
					warning('%s not found in %s' % (samp, opt.pops))
					prefixes = [os.path.commonprefix([s, samp]) for s in samppop.keys()]
					longpref = max(prefixes, key=len)
					if len(longpref) > 6: #ARBITRARY HACK TO FIND BEST MATCHING SAMPLE
						match = samppop.keys()[prefixes.index(longpref)]
						samppop[samp] = samppop[match]
						warning('Using %s\t%s' % (match, samppop[match]))

				if not samppop[samp] in poplist: # error if samp not in opt.pops
					poplist.append(samppop[samp])
					popsampnums[samppop[samp]] = []
				popsampnums[samppop[samp]].append(sampn[samp])
		else:
			poplist = ['all_samples']
			popsampnums['all_samples'] = range(nsamp)

		if opt.condsamp:
			for i, name in enumerate(sampnames):
				if name.find(condsamp) >= 0:
					condcol = i + VCF_FIXEDCOLS
#				   print(i, condsamp, tok[condcol])
					break
		break


if opt.rates:
	if opt.depths:
		depcount = [0 for i in range(nsamp)]
	hetcount = [0 for i in range(nsamp)]
	hnrcount = [0 for i in range(nsamp)]
	nsites = 0
	firstcoord = ''

if opt.psmcfa:
	faout = FastaStream(sys.stdout)

#if opt.header:
#	if opt.pops:
#		fout.write('#VCF_POPULATIONS:\t%s\n' % ('\t'.join(poplist)))
#	if not opt.no_calls:
#		fout.write('\t'.join(['#SAMPLES'] + sampnames) + '\n')
if opt.header: #column headers
	header = []
	header += ['CHR', 'POS']
	if opt.depthsonly:
		header += sampnames
	if opt.aims:
		header += ['***HEADER_INFO_TODO***']
	if not opt.no_calls:
		header += sampnames
	if opt.segsites:
		header += ['SEGSITE']
	if opt.diversity:
		for pop in poplist:
			header += [pop + '_N', pop + '_AC', pop + '_AF']
	if opt.Fst:
		header += ['Fst']
	fout.write('\t'.join(header) + '\n')


curchr = ''
for line in fin:
	if line.startswith('#'):
#		sys.stdout.write(line)
		continue

	tok = line.split()

	if tok[0] != curchr:
		sep = 0
		lastpos = 0
		curchr = tok[0]
		if opt.psmcfa:
			nextbin = PSMC_BINSIZE
			faout = FastaStream(sys.stdout)
			faout.newrec(curchr)
			inmask = False
			nN = 0
			outchar = 'T'

	if opt.replacecalls and tok[0]+tok[1] in reprecord:
		tok[VCF_FIXEDCOLS-1:] = reprecord[tok[0]+tok[1]]
	
	if opt.vars and tok[4] == '.':
		continue

#	sys.stdout.write(line)
	if not opt.indels and (len(tok[3]) > 1 or re.search('INDEL', tok[7])):
		continue

#	if len(tok[4]) > 1: # multiallelic sites
#		continue

	if opt.condsamp and not tok[condcol].startswith('0/1'): # condition on variant in cond_sample
		continue

	if opt.depths or opt.mediandep:
		formats = tok[VCF_FIXEDCOLS - 1].split(':')
		depfield = formats.index('DP')
		dep = [x.split(':')[depfield] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]

	if opt.depthsonly:
		fout.write('\t'.join(tok[0:2] + dep) + '\n')
		continue

	outvals = tok[0:2]

	if opt.haploid:
		if (tok[4] == "."):
			varvals = ['0' for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
			segsite = False
		else:
			varvals = [x[0] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
		vargts = varvals
	else:
		if (tok[4] == "."):
			varvals = ['00' for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
			vargts = ['0/0' for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
			segsite = False
		else:
			varvals = [x[0] + x[2] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
			vargts = [x[0:3] for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]

		if opt.pseudodip:
			hetsite = False
			for gt in vargts:
				if gt[0] != gt[2]:
					hetsite = True
					break
			if hetsite:
				continue

	if opt.pseudodip:
# combine consecutive varvals elements into new varvals (assumes none are het)
# TODO?: only check/exclude hets in consecutive pairs
		varvals = [varvals[i][0] + varvals[i+1][0] for i in range(0, len(varvals) - 1, 2)]
		vargts = [vargts[i][0] + '|' + vargts[i+1][0] for i in range(0, len(vargts) - 1, 2)]
#		debug(vargts)

# determine if segregating site
# ideally use FQ > 0 but vcf-subset doesn't currently adjust this
#		fq = re.search('FQ=(-?[0-9.]+)', tok[7])
#		if fq > 0: #het
#		sys.stdout.write(line)
	varvalstr = '-'.join(varvals)
#	debug(valstr)
	segsite = '0' in varvalstr and '1' in varvalstr

	if opt.alleles: # convert 0, 1 to alleles in output
		sitealleles = [tok[3]] + tok[4].split(',')
		trtab = maketrans(''.join([str(x) for x in range(len(sitealleles))]), ''.join(sitealleles))
		varalleles = [x.translate(trtab) for x in varvals]  
		varvals = varalleles  
		varvalstr = '-'.join(varvals)
		vargts = [x.translate(trtab) for x in vargts]  
#		debug([tok[1], sitealleles, vargts])

	if opt.aims:
		valpops = {}
		popvals = {}
		for isn, samp in enumerate(sampnames):
			spop = samppop[samp]
			for val in varvals[isn]:
				if val in valpops:
					valpops[val].add(spop)
				else:
					valpops[val] = set([spop])
				if spop in popvals:
					popvals[spop].add(val)
				else:
					popvals[spop] = set([val])
		for val in valpops:
			if len(valpops[val]) == 1: # require uniqueness of marker to pop
				spop = valpops[val].pop() 
				if len(popvals[spop]) == 1: # require fixation of marker in pop
					debug(varvals)
#					debug(valpops)
					debug(popvals)
					spvals = ['/'.join(list(popvals[x])) for x in poplist]
					debug(spvals)
					fout.write('\t'.join(outvals + spvals + [spop]) + '\n')
#					fout.write('\t'.join(outvals + [val] + spvals + [spop]) + '\n')
		continue

	if opt.segsep:
		pos = int(tok[1])

		while pos > lastpos:
			if opt.callmask and lastpos in callmask:
#				print('pos = %d; Masking %d to %d' % (pos, lastpos, callmask[lastpos]))
				lastpos = callmask[lastpos]
			else:
				sep += 1
			lastpos += 1
		if opt.callmask and lastpos in callmask:
#			print('pos = %d; Masking %d to %d' % (pos, lastpos, callmask[lastpos]))
			lastpos = callmask[lastpos] + 1
		elif lastpos == pos and segsite:
			if len(vargts) > 1:
				allstr = ','.join(phasecombs(vargts))
				fout.write('\t'.join(outvals + [str(sep)] + [allstr]) + '\n')
			else:
				fout.write('\t'.join(outvals + [str(sep)] + varvals) + '\n')
			sep = 0
		continue

	elif opt.psmcfa:
		if segsite:
			pos = int(tok[1])

#			print('site at %d' % pos)
			while lastpos < nextbin and pos > nextbin - PSMC_BINSIZE:
				lastpos += 1
				if opt.callmask and lastpos in callmask:
#					print('Mask %d to %d : %d' % (lastpos, callmask[lastpos], callmask[lastpos]-lastpos + 1))
					inmask = True
					endmask = callmask[lastpos]
				if inmask:
					nN += 1
					if lastpos == endmask:
						inmask = False
				if lastpos == pos:
					outchar = 'K'
				if lastpos == nextbin:
					if float(nN) / PSMC_BINSIZE > PSMC_N_RATIO:
						outchar = 'N'
					faout.write(outchar)
#					print(' %d-%d: %d Ns, pos = %d' % (nextbin - PSMC_BINSIZE + 1, nextbin, nN, pos))
					nextbin += PSMC_BINSIZE
					outchar = 'T'
					nN = 0
		continue

	if opt.rates:
		if opt.depths:
			depcount = [depcount[i] + int(dep[i]) for i in range(nsamps)]
		hetcount = [hetcount[i] + (varvals[i] == '01') for i in range(nsamps)]
		hnrcount = [hnrcount[i] + (varvals[i] == '11') for i in range(nsamps)]
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

		alcount_total = 0
		samps_total = 0
		Hs_mean = 0
		for pop in poplist:
			pvarvalstr = '-'.join([varvals[i] for i in popsampnums[pop]])
			if not opt.no_calls:
				outvals.append(pvarvalstr)

			if opt.diversity or opt.Fst:
	#			alfreq = eval(re.search('AF1=([0-9.e\-]+)', tok[7]).group(1))# eval since vcf uses sci notation
	#			alcount = int(re.search('AC1=([0-9]+)', tok[7]).group(1))
	#			neidiv = 1 - alfreq**2 - (1 - alfreq)**2
	#			outvals += [str(nsamp), str(alfreq), str(neidiv), str(alcount), str(alcount/float(nsamp))]
				alcount = pvarvalstr.count('1')
				pnsamp = len(popsampnums[pop])
				pfreq = alcount/float(2 * pnsamp)
				alcount_total += alcount
				samps_total += pnsamp
				Hs_mean += 2 * pfreq * (1 - pfreq)
				if opt.diversity:
					outvals += [str(pnsamp), str(alcount), str(pfreq)]

		if opt.Fst:
			Hs_mean = Hs_mean / len(poplist)
			tfreq = alcount_total/float(2 * samps_total)
			Ht = 2 * tfreq * (1 - tfreq)
			Fst = (Ht - Hs_mean) / Ht
#			outvals += [str(Hs_mean), str(Ht), str(Fst)]
			outvals += [str(Fst)]

		fout.write('\t'.join(outvals) + '\n')

if opt.rates:
	fout.write('\t'.join(['#REGION'] + [firstcoord + '-' + lastcoord]) + '\n')
	fout.write('\t'.join(['#NSITES'] + [str(nsites)]) + '\n')
	if opt.depths:
		fout.write('\t'.join(['#MEANDEP'] + [str(float(x)/nsites) for x in depcount]) + '\n')
	fout.write('\t'.join(['#HETCOUNT'] + [str(x) for x in hetcount]) + '\n')
	fout.write('\t'.join(['#HETRATE'] + [str(float(x)/nsites) for x in hetcount]) + '\n')
	fout.write('\t'.join(['#HNRCOUNT'] + [str(x) for x in hnrcount]) + '\n')
	fout.write('\t'.join(['#HNRRATE'] + [str(float(x)/nsites) for x in hnrcount]) + '\n')
	fout.write('\t'.join(['#HNR/HET'] + [str(float(hnrcount[i])/hetcount[i]) for i in range(nsamp)]) + '\n')
#	fout.write('\t'.join(outvals) + '\n')
