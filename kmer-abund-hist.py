#! /usr/bin/env python
"""
Plot subsampled k-mer abundance histograms from sourmash sigs.

See https://github.com/dib-lab/sourmash/issues/910 for sourmash integration.

Relies on sourmash >= 4.0.0rc1, termplotlib
"""
import sys
import argparse
import numpy
import collections
import termplotlib as tpl

import sourmash
from sourmash import sourmash_args
from sourmash.cli.utils import add_moltype_args, add_ksize_arg
from sourmash.logging import set_quiet, notify


def main():
    p = argparse.ArgumentParser()
    p.add_argument('signatures', nargs='+')
    p.add_argument(
        '-q', '--quiet', action='store_true',
        help='suppress non-error output'
    )
    p.add_argument(
        '-o', '--output', metavar='FILE',
        help='output histogram to this file (in CSV format)'
    )
    p.add_argument(
        '--abundances', metavar='FILE',
        help='output hashes and abundances to this file (in CSV format)')
    p.add_argument(
        '--md5', default=None,
        help='select signatures whose md5 contains this substring'
    )
    p.add_argument(
        '--name', default=None,
        help='select signatures whose name contains this substring'
    )
    p.add_argument(
        '--max', type=int, default=None,
        help='max value for histogram range (default none)')
    p.add_argument(
        '--min', type=int, default=None,
        help='min value for histogram range (default none)')
    p.add_argument(
        '--bins', type=int, default=10,
        help='number of bins (default 10)')
    add_ksize_arg(p, 31)
    add_moltype_args(p)

    args = p.parse_args()
    return abundhist(args)


def abundhist(args):
    """
    output abundance histogram and/or raw abundances.
    """

    set_quiet(args.quiet)
    moltype = sourmash_args.calculate_moltype(args)

    outlist = []
    total_loaded = 0
    for filename in args.signatures:
        siglist = sourmash.load_file_as_signatures(filename, ksize=args.ksize,
                                                   select_moltype=moltype)
        siglist = list(siglist)

        total_loaded += len(siglist)

        # select!
        if args.md5 is not None:
            siglist = [ ss for ss in siglist if args.md5 in ss.md5sum() ]
        if args.name is not None:
            siglist = [ ss for ss in siglist if args.name in ss.name() ]

    notify("loaded {} total that matched ksize & molecule type",
           total_loaded)
    if len(siglist) != total_loaded:
        notify("selected {} via name / md5 selectors".format(len(siglist)))
    notify('')

    counts_d = collections.defaultdict(int)
    for ss in siglist:
        for hashval, abund in ss.minhash.hashes.items():
            counts_d[hashval] += abund

    all_counts = list(counts_d.values())

    min_range = 1
    if args.min is not None:
        min_range = args.min
    max_range = max(all_counts)
    if args.max is not None:
        max_range = args.max

    n_bins = args.bins
    if max_range - min_range + 1 < n_bins:
        n_bins = max_range - min_range + 1

    # make hist
    counts, bin_edges = numpy.histogram(all_counts,
                                        range=(min_range, max_range),
                                        bins=n_bins)
    bin_edges = bin_edges.astype(int)

    # plot
    fig = tpl.figure()
    f = fig.barh(counts, [ str(x) for x in bin_edges[1:] ], force_ascii=True)
    fig.show()

    # output histogram in csv?
    if args.output:
        with FileOutput(args.output, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['count', 'n_count'])
            for nc, c in zip(counts, bin_edges[1:]):
                w.writerow([c, nc])

    # output raw counts tagged with hashval?
    if args.abundances:
        with FileOutput(args.abundances, 'wt') as fp:
            w = csv.writer(fp)
            w.writerow(['hashval', 'count'])
            for hashval, count in counts_d.items():
                w.writerow([hashval, count])


if __name__ == '__main__':
    main()
