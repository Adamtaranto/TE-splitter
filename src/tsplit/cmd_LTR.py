from tsplit.logs import init_logging
from tsplit.utils import segWrite, tSplitchecks, check_tools
from tsplit.parseAlign import getLTRs

def parse_ltr_args(args):
    return args

def main(args=None):
    """Do the work."""
    if args is None:
        args = parse_ltr_args()

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

    # If LTR mode, search for terminal repeats on same strand
    segments = getLTRs(
        args.infile,
        flankdist=args.maxdist,
        minterm=args.minterm,
        minseed=args.minseed,
        minid=args.minid,
        diagfactor=args.diagfactor,
        report=args.splitmode,
        temp=args.outdir,
        alignTool=args.method,
        keeptemp=args.keeptemp,
    )

    segWrite(outpath, segs=segments)

if __name__ == "__main__":
    main()