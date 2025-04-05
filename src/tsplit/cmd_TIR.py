from tsplit.logs import init_logging
from tsplit.parseAlign import getTIRs
from tsplit.utils import check_tools, segWrite, tSplitchecks


def main(args=None):
    """Do the work."""
    # Set up logging
    init_logging(loglevel=args.loglevel)

    # Check for required programs.
    required_tools = ['delta-filter', 'nucmer', 'show-coords']
    optional_tools = ['blastn']

    if args.method == 'blastn':
        required_tools.append('blastn')
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
