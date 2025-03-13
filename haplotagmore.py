#!/usr/bin/env python

import sys
import argparse
import pysam
import numpy as np

from collections import defaultdict as dd
from cyvcf2 import VCF
from bx.intervals.intersection import Intersecter, Interval

from tqdm import tqdm

import logging
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_het_vars(args):
    vcf = VCF(args.vcf)

    p_idx = 0

    snps = dd(Intersecter)
    count = 0

    for rec in tqdm(vcf):

        if rec.FILTER is not None:
            if not args.all_snps:
                continue
        
        if len(rec.ALT) != 1:
            continue

        ph1 = False
        ph2 = False

        p_gt = rec.genotypes[p_idx][:2]

        if ph2 == ph1 == False:
            if rec.gt_phases[p_idx]: # was previously phased
                if 0 in p_gt and 1 in p_gt:
                    if p_gt[0] == 1:
                        ph1 = True
                    
                    if p_gt[1] == 1:
                        ph2 = True

        category = None

        assert not ph1 == ph2 == True

        if ph1 != ph2:
            category = 'ph2'
            if ph1:
                category = 'ph1'

        if category is None:
            continue

        snps[rec.CHROM].add_interval(Interval(rec.POS, rec.POS, value=(rec.REF, rec.ALT[0], category)))
        count += 1
    
    logger.info('found %d heterozygous variants inherited from one parent' % count)

    return snps


def get_indels(read, genome):
    indels = {}

    query_lookup = dict([(q, r) for q, r in read.get_aligned_pairs() if q is not None and r is not None])

    q_pos = 0
    lastop = None
    lastref = None
    for op in read.cigartuples:
        assert op[0] < 6

        if op[0] == 4: # soft clip
            q_pos += op[1]
        
        if op[0] == 0: # match
            q_pos += op[1]

        if op[0] == 3: # skip
            pass

        if op[0] == 5: # hardclip
            pass
        
        if op[0] == 1: # ins
            if q_pos-1 not in query_lookup:
                logger.warning('position %d not in query_lookup read id %s' % (q_pos-1, read.qname))
                continue

            ref_start = query_lookup[q_pos-1]
            lastref = ref_start
            bases = read.seq[q_pos:q_pos+op[1]]
            indels[ref_start] = ('I', bases)
            q_pos += op[1]

        
        if op[0] == 2: # del
            # try:
            if lastop == 1:
                ref_start = lastref
            else:
                try:
                    ref_start = query_lookup[q_pos-1]
                except KeyError:
                    ref_start = query_lookup[q_pos]
                finally:
                    return indels
                    
            lastref = ref_start
            if read.reference_name in genome.references:
                bases = genome.fetch(read.reference_name, ref_start, ref_start+op[1])
                indels[ref_start] = ('D', bases)
            else:
                logger.warning(f'{read.reference_name} not in genome')

        lastop = op[0]

    return indels


def main(args):
    if args.debug:
        logger.setLevel(logging.DEBUG)

    logger.info('%s fetch informative variants from %s...' % (args.bam, args.vcf))
    snps = get_het_vars(args)

    if len(snps) == 0:
        sys.exit('no phased het snps in vcf, exiting.')

    bam = pysam.AlignmentFile(args.bam)
    genome = pysam.Fastafile(args.ref)

    out_fn = '.'.join(args.bam.split('.')[:-1]) + '.tagmore.bam'

    logger.info('output bam: %s' % out_fn)
    out = pysam.AlignmentFile(out_fn, 'wb', template=bam)

    stats_votes = []
    stats_margin = []
    count_reads = 0
    count_success = 0
    already_tagged = 0

    votes = {}

    logger.info('assign variants to reads...')

    for read in tqdm(bam.fetch(), total=bam.mapped+bam.unmapped):
        if read.is_duplicate:
            continue

        if not read.seq:
            continue

        if read.is_unmapped:
            out.write(read)
            continue

        count_reads += 1

        # already phased?

        if read.has_tag('HP'):
            already_tagged += 1
            continue
        
        elif read.has_tag('PS'):
            raise ValueError(f'malformed record ({read.query_name}) has PS tag without HP tag')
        
        base_lookup = dict([(r, read.seq[q]) for q, r in read.get_aligned_pairs() if q is not None])

        indels = get_indels(read, genome)

        if read.qname not in votes:
            votes[read.qname] = [0,0]

        matched = 0
        unmatched = 0
        ins_matched = 0
        del_matched = 0

        for snp in snps[read.reference_name].find(read.reference_start, read.reference_end):
            ref, alt, origin = snp.value
            if len(ref) == len(alt) == 1: #snp
                if snp.start-1 in base_lookup:
                    b = base_lookup[snp.start-1]

                    if b.upper() == alt.upper():
                        if origin == 'ph2':
                            votes[read.qname][1] += 1
                        if origin == 'ph1':
                            votes[read.qname][0] += 1
                        matched += 1

                    elif b.upper() == ref.upper(): # inverse of alt
                        if origin == 'ph2':
                            votes[read.qname][0] += 1
                        if origin == 'ph1':
                            votes[read.qname][1] += 1
                        matched += 1
                    
                    else:
                        unmatched += 1

            if args.ignore_indels:
                continue

            elif len(ref) == 1: # ins
                if snp.start-1 in indels:
                    if indels[snp.start-1][0] == 'I':
                        if args.loose_indels:
                            min_l = min(len(alt[1:]), len(indels[snp.start-1][1]))
                            if alt[1:][-min_l:] == indels[snp.start-1][1][-min_l:]:
                                ins_matched += 1
                                if origin == 'ph2':
                                    votes[read.qname][1] += 1
                                if origin == 'ph1':
                                    votes[read.qname][0] += 1

                        else:
                            if alt[1:] == indels[snp.start-1][1]:
                                ins_matched += 1
                                if origin == 'ph2':
                                    votes[read.qname][1] += 1
                                if origin == 'ph1':
                                    votes[read.qname][0] += 1

            elif len(alt) == 1: # del
                if snp.start-1 in indels:
                    if indels[snp.start-1][0] == 'D':
                        if args.loose_indels:
                            min_l = min(len(ref), len(indels[snp.start-1][1]))
                            if ref[:min_l] == indels[snp.start-1][1][:min_l]:
                                del_matched += 1
                                if origin == 'ph2':
                                    votes[read.qname][1] += 1
                                if origin == 'ph1':
                                    votes[read.qname][0] += 1
                        
                        else:
                            if ref == indels[snp.start-1][1]:
                                del_matched += 1
                                if origin == 'ph2':
                                    votes[read.qname][1] += 1
                                if origin == 'ph1':
                                    votes[read.qname][0] += 1

    bam.reset()
    
    logger.info('tagging reads...')

    last_PS = None

    for read in tqdm(bam.fetch(), total=bam.mapped+bam.unmapped):
        if read.qname not in votes:
            if not args.phased_only:
                out.write(read)
            continue

        if read.has_tag('PS'):
            last_PS = read.get_tag('PS')

        minvotes = int(args.minvotes)
        minmargin = int(args.minmargin)
     
        if sum(votes[read.qname]) < minvotes:
            logger.debug('less than %d votes for phased read %s' % (minvotes, read.qname))
            if not args.phased_only:
                out.write(read)
            continue

        if abs(votes[read.qname][0]-votes[read.qname][1]) < minmargin:
            logger.debug('less than %d vote margin for phased read %s' % (minmargin, read.qname))
            if not args.phased_only:
                out.write(read)
            continue

        if votes[read.qname][0] == votes[read.qname][1]:
            logger.debug('ambiguous origin for read %s' % read.qname)
            if not args.phased_only:
                out.write(read)
            continue

        hp = None

        if votes[read.qname][0] > votes[read.qname][1]: # ph1
            hp = 1

        if votes[read.qname][1] > votes[read.qname][0]: # ph2
            hp = 2

        stats_votes.append(sum(votes[read.qname]))
        stats_margin.append(abs(votes[read.qname][0]-votes[read.qname][1]))

        count_success += 1

        logger.debug('read %s (%s) matched: %d, unmatched: %d, ins_matched: %d, del_matched: %d, vote_1: %d, vote_2: %d, hp: %d' % 
                    (read.qname, str(read.is_reverse), matched, unmatched, ins_matched, del_matched, votes[read.qname][0], votes[read.qname][1], hp))

        if last_PS is None:
            last_PS = read.reference_start

        read.set_tag('PS', last_PS)
        read.set_tag('HP', hp)
        read.set_tag('xH', 'HTM')
        out.write(read)
    
    out.close()

    logger.info('%s: tagged before haplotagmore: %d' % (args.bam, already_tagged))
    logger.info('%s: tagged after haplotagmore: %d' % (args.bam, already_tagged+count_success))
    logger.info('%s: tagged %d alignments out of %d untagged alignments' % (args.bam, count_success, count_reads-already_tagged))
    logger.info('%s: mean votes: %.2f, median votes: %.2f' % (args.bam, np.mean(stats_votes), np.median(stats_votes)))
    logger.info('%s: mean margin: %.2f, median margin: %.2f' % (args.bam, np.mean(stats_margin), np.median(stats_margin)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='assign alignments to parent of origin based on snps and indels')
    parser.add_argument('-v', '--vcf', required=True, help='vcf')
    parser.add_argument('-b', '--bam', required=True, help='bam')
    parser.add_argument('-r', '--ref', required=True, help='indexed fasta')
    parser.add_argument('--minvotes', default=1, help='default = 1')
    parser.add_argument('--minmargin', default=1, help='default = 1')
    parser.add_argument('--debug', action='store_true', default=False)
    parser.add_argument('--all_snps', action='store_true', default=False, help='ignore vcf filters (default: PASS only)')
    parser.add_argument('--loose_indels', action='store_true', default=False, help='allow left end of deletion or right end of insertion matching (useful for nanopore indels)')
    parser.add_argument('--ignore_indels', action='store_true', default=False, help='do not consider indels')
    parser.add_argument('--phased_only', action='store_true', default=False, help='only output phased reads')
    args = parser.parse_args()
    main(args)
