import argparse
import sys
from tsplit.cmd_TIR import main as tir_main
from tsplit.cmd_LTR import main as ltr_main


def main():
    parser = argparse.ArgumentParser(
        description="Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs)."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # TIR subcommand
    tir_parser = subparsers.add_parser(
        "TIR", help="Extract terminal repeats from DNA transposons (TIRs)."
    )
    tir_parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'INFO'",
    )
    tir_parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        default=None,
        help="Multifasta containing complete elements.",
    )
    tir_parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="All output files begin with this string. (Default: [infile name])",
    )
    tir_parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        default=None,
        help="Write output files to this directory. (Default: cwd)",
    )
    tir_parser.add_argument(
        "--splitmode",
        default="split",
        choices=["all", "split", "internal", "external", None],
        help='all= Report input sequence as well as internal and external segments. split= Report internal and external segments after splitting. internal = Report only internal segments external = Report only terminal repeat segments. If set to "None" then only synthetic MITES will be reported if --makemites is also set. (Default: split)',
    )
    tir_parser.add_argument(
        "--makemites",
        action="store_true",
        default=False,
        help="Attempt to construct synthetic MITE sequences from TIRs.",
    )
    tir_parser.add_argument(
        "--keeptemp",
        action="store_true",
        default=False,
        help="If set do not remove temp directory on completion.",
    )
    tir_parser.add_argument(
        "-m",
        "--maxdist",
        type=int,
        default=10,
        help="Terminal repeat candidates must be no more than this many bases from end of input element. (Default: 10) Note: Increase this value if you suspect that your element is nested within some flanking sequence.",
    )
    tir_parser.add_argument(
        "--minid",
        type=float,
        default=80.0,
        help="Minimum percentage identity between terminal repeat pairs. As float. (Default: 80.0)",
    )
    tir_parser.add_argument(
        "--minterm",
        type=int,
        default=10,
        help='Minimum length for a terminal repeat to be considered. Equivalent to nucmer "--mincluster" (Default: 10)',
    )
    tir_parser.add_argument(
        "--minseed",
        type=int,
        default=5,
        help='Minimum length of a maximal exact match to be included in final match cluster. Equivalent to nucmer "--minmatch". (Default: 5)',
    )
    tir_parser.add_argument(
        "--diagfactor",
        type=float,
        default=0.2,
        help="Maximum diagonal difference factor for clustering of matches within nucmer, i.e. diagonal difference / match separation (default 0.20) Note: Increase value for greater tolerance of indels between terminal repeats.",
    )
    tir_parser.add_argument(
        "--method",
        default="blastn",
        choices=["blastn", "nucmer"],
        help='Select alignment method: "blastn" or "nucmer".(Default: blastn)',
    )

    # LTR subcommand
    ltr_parser = subparsers.add_parser(
        "LTR", help="Extract terminal repeats from retrotransposons (LTRs)."
    )
    ltr_parser.add_argument(
        "--loglevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level. Default: 'DEBUG'",
    )
    ltr_parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        default=None,
        help="Multifasta containing complete elements.",
    )
    ltr_parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        default=None,
        help="All output files begin with this string. (Default: [infile name])",
    )
    ltr_parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        default=None,
        help="Write output files to this directory. (Default: cwd)",
    )
    ltr_parser.add_argument(
        "--splitmode",
        default="split",
        choices=["all", "split", "internal", "external", None],
        help='all= Report input sequence as well as internal and external segments. split= Report internal and external segments after splitting. internal = Report only internal segments external = Report only terminal repeat segments. If set to "None" then only synthetic MITES will be reported if --makemites is also set. (Default: split)',
    )
    ltr_parser.add_argument(
        "--keeptemp",
        action="store_true",
        default=False,
        help="If set do not remove temp directory on completion.",
    )
    ltr_parser.add_argument(
        "-m",
        "--maxdist",
        type=int,
        default=10,
        help="Terminal repeat candidates must be no more than this many bases from end of input element. (Default: 10) Note: Increase this value if you suspect that your element is nested within some flanking sequence.",
    )
    ltr_parser.add_argument(
        "--minid",
        type=float,
        default=80.0,
        help="Minimum percentage identity between terminal repeat pairs. As float. (Default: 80.0)",
    )
    ltr_parser.add_argument(
        "--minterm",
        type=int,
        default=10,
        help='Minimum length for a terminal repeat to be considered. Equivalent to nucmer "--mincluster" (Default: 10)',
    )
    ltr_parser.add_argument(
        "--minseed",
        type=int,
        default=5,
        help='Minimum length of a maximal exact match to be included in final match cluster. Equivalent to nucmer "--minmatch". (Default: 5)',
    )
    ltr_parser.add_argument(
        "--diagfactor",
        type=float,
        default=0.2,
        help="Maximum diagonal difference factor for clustering of matches within nucmer, i.e. diagonal difference / match separation (default 0.20) Note: Increase value for greater tolerance of indels between terminal repeats.",
    )
    ltr_parser.add_argument(
        "--method",
        default="blastn",
        choices=["blastn", "nucmer"],
        help='Select alignment method: "blastn" or "nucmer".(Default: blastn)',
    )

    args = parser.parse_args()

    if args.command == "TIR":
        tir_main(args)
    elif args.command == "LTR":
        ltr_main(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
