import argparse

from tsplit.logs import init_logging
from tsplit.utils import segWrite, tSplitchecks, check_tools
from tsplit.parseAlign import getTIRs


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Extract terminal repeats from DNA transposons (TIRs). \
                    Optionally, compose synthetic MITES from complete DNA transposons.",
        prog="tsplit-TIR",
    )
    parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'DEBUG'",
    )
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        default=None,
        help="Multifasta containing complete elements.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="All output files begin with this string. \
                                        (Default: [infile name])",
    )
    parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        default=None,
        help="Write output files to this directory. \
                                        (Default: cwd)",
    )
    parser.add_argument(
        "--splitmode",
        default="split",
        choices=["all", "split", "internal", "external", None],
        help='all= Report input sequence as well as internal and external segments. \
                                        split= Report internal and external segments after splitting. \
                                        internal = Report only internal segments \
                                        external = Report only terminal repeat segments. \
                                        If set to "None" then only synthetic MITES will be reported if --makemites is also set. \
                                        (Default: split)',
    )
    parser.add_argument(
        "--makemites",
        action="store_true",
        default=False,
        help="Attempt to construct synthetic MITE sequences from TIRs.",
    )
    parser.add_argument(
        "--keeptemp",
        action="store_true",
        default=False,
        help="If set do not remove temp directory on completion.",
    )
    parser.add_argument(
        "-m",
        "--maxdist",
        type=int,
        default=10,
        help="Terminal repeat candidates must be no more than this many bases from end of input element. (Default: 10)\
                                        Note: Increase this value if you suspect that your element is nested within some flanking sequence.",
    )
    parser.add_argument(
        "--minid",
        type=float,
        default=80.0,
        help="Minimum percentage identity between terminal repeat pairs. As float. \
                                        (Default: 80.0)",
    )
    parser.add_argument(
        "--minterm",
        type=int,
        default=10,
        help='Minimum length for a terminal repeat to be considered. \
                                        Equivalent to nucmer "--mincluster" \
                                        (Default: 10)',
    )
    parser.add_argument(
        "--minseed",
        type=int,
        default=5,
        help='Minimum length of a maximal exact match to be included in final match cluster. \
                                        Equivalent to nucmer "--minmatch". \
                                        (Default: 5)',
    )
    parser.add_argument(
        "--diagfactor",
        type=float,
        default=0.2,
        help="Maximum diagonal difference factor for clustering of matches within nucmer, \
                                        i.e. diagonal difference / match separation (default 0.20) \
                                        Note: Increase value for greater tolerance of indels between terminal repeats.",
    )
    parser.add_argument(
        "--method",
        default="blastn",
        choices=["blastn", "nucmer"],
        help='Select alignment method: "blastn" or "nucmer".(Default: blastn)',
    )
    args = parser.parse_args()
    return args


def main():
    """Do the work."""

    # Get cmd line args
    args = mainArgs()

    # Set up logging
    init_logging(loglevel=args.loglevel)

    # Check for required programs.
    required_tools = ["delta-filter", "nucmer", "show-coords"]
    optional_tools = ["blastn"]

    if args.method == "blastn":
        required_tools.append("blastn")
        optional_tools = []

    check_tools(required_tools=required_tools, optional_tools=optional_tools)

    # Create output paths as required
    outpath = tSplitchecks(args)

    # Search for inverted terminal repeats
    # Optionally construct synthetic MITE from TIRs if detected
    segments = getTIRs(
        args.infile,
        flankdist=args.maxdist,
        minterm=args.minterm,
        minseed=args.minseed,
        minid=args.minid,
        diagfactor=args.diagfactor,
        mites=args.makemites,
        report=args.splitmode,
        alignTool=args.method,
        temp=args.outdir,
        keeptemp=args.keeptemp,
    )

    segWrite(outpath, segs=segments)


if __name__ == "__main__":
    main()
