#!/software/bin/python
# Aylwyn Scally 2013

import sys
import optparse
#import os
import re
#import os.path
#from numpy import median
import logging
from logging import error, warning, info, debug, critical

loglevel = logging.WARNING

# defaults
VCF_FIXEDCOLS = 9
#HWEMIN = 0.05

p = optparse.OptionParser()
p.add_option('-H', '--noheader', action='store_true', default = False)
p.add_option('-d', '--mindep', dest = 'mindeparg', default = '10', help = 'min per-individual raw read depth [%default]')
p.add_option('-D', '--maxdep', dest = 'maxdeparg', default = '100', help = 'max per-individual raw read depth [%default]')
p.add_option('-n', '--min_passed', dest = 'minpass', type = 'float', default = -1.0, help = 'min samples passing depth filters: n < 0: all samples; 0 < n < 1: fraction of samples; n >= 1: number of samples')
p.add_option('-q', '--Qmin', type = 'float', default = 20.0, help = 'min consensus quality [%default]')
p.add_option('--MQmin', default = 10, help = 'min r.m.s. mapping quality of covering reads [%default]')
p.add_option('--MQ0maxfrac', type = 'float', default = 0.25, help = 'max raction of reads with zero mapping quality [%default]')
#p.add_option('--FQmin', type = 'float', default = 10.0, help = '')
p.add_option('--PV4min', type = 'float', default = 0.0001, help = 'min p-value for strand bias, baseQ bias, mapQ bias and tail distance bias [%default]')
p.add_option('--HWEmin', type = 'float', default = 0.0, help = 'min HWE chi^2 test p-value [%default]')
#DP4: number of alternate bases > 2


opt, args = p.parse_args()
mindep = [int(x) for x in opt.mindeparg.split(',')]
maxdep = [int(x) for x in opt.maxdeparg.split(',')]

if args:
	fin = open(args[0])
else:
	fin = sys.stdin
fout = sys.stdout

for line in fin:
	if not opt.noheader:
		sys.stdout.write(line)
	if line.startswith('#CHROM'):
		tok = line.split()
		nsamp = len(tok) - VCF_FIXEDCOLS
#		sampnames = tok[VCF_FIXEDCOLS:]
		break

if opt.minpass < 0:
	opt.minpass = nsamp
elif opt.minpass > 0 and opt.minpass < 1:
	opt.minpass = int(opt.minpass * nsamp)
elif opt.minpass >= 1 and opt.minpass <= nsamp:
	opt.minpass = int(opt.minpass)
elif opt.minpass > nsamp:
	error('%d samples, %d given for min number passing depth filter' % (nsamp, opt.minpass))
	sys.exit(2)

if len(mindep) > 1 and nsamp > len(mindep):
	error('%d samples, %d values given for min depth' % (nsamp, len(mindep)))
	sys.exit(2)
elif len(mindep) == 1:
	mindep += [mindep[0]] * (nsamp - 1)

if len(maxdep) > 1 and nsamp > len(maxdep):
	error('%d samples, %d values given for max depth' % (nsamp, len(maxdep)))
	sys.exit(2)
elif len(maxdep) == 1:
	maxdep += [maxdep[0]] * (nsamp - 1)

#depfield = 2
for line in fin:
	if line.startswith('#'):
		sys.stdout.write(line)
		continue
	tok = line.split()

#	if len(tok[4]) > 1: # multiallelic sites
#		continue

	if float(tok[5]) < opt.Qmin:
		continue

	formats = tok[VCF_FIXEDCOLS - 1].split(':')

	depfield = formats.index('DP')
	dp = [int(x.split(':')[depfield]) for x in tok[VCF_FIXEDCOLS:(VCF_FIXEDCOLS+nsamp)]]
	passed = sum([dp[i] >= mindep[i] and dp[i] <= maxdep[i] for i in range(nsamp)])
	if passed < opt.minpass:
		continue

	if int(re.search('MQ=([0-9]+)', tok[7]).group(1)) < opt.MQmin:
		continue

	DP = int(re.search('DP=([0-9]+)', tok[7]).group(1))
	if int(re.search('MQ0=([0-9]+)', tok[7]).group(1)) > opt.MQ0maxfrac * DP:
		continue

#	if abs(float(re.search('FQ=(-?[0-9.]+)', tok[7]).group(1))) < opt.FQmin:
#		continue

	filt = False
	pv4p = re.search('PV4=([0-9.e\-,]+)', tok[7])
	if pv4p:
		for pv in pv4p.group(1).split(','):
			if eval(pv) < opt.PV4min:
				filt = True
	if filt:
		continue

#	DP4: number of alternate bases > 2 ONLY APPLIES TO VARIANTS

	if nsamp > 1:
		if opt.HWEmin > 0.0:
			hwep = re.search('HWE=([0-9.e\-]+)', tok[7])
			if hwep and eval(hwep.group(1)) < opt.HWEmin:
				continue

	sys.stdout.write(line)

#	alfreq = eval(re.search('AF1=([0-9.e\-]+)', tok[7]).group(1))# eval since vcf uses sci notation
#	alcount = int(re.search('AC1=([0-9]+)', tok[7]).group(1))
