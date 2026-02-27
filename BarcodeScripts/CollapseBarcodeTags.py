import difflib
import re
import sys
from argparse import ArgumentParser, Namespace
from collections import defaultdict
from datetime import timedelta
from subprocess import PIPE, Popen, CalledProcessError, run
from time import time as now

from Levenshtein import hamming, distance as levenshtein, opcodes

# Regular expression for detecting raw barcodes in bam file
RG_RE = re.compile(r"RG:Z:([^\s]+)")


def main(args):
    def elapsed():
        nonlocal start
        return timedelta(seconds=int(now() - start))

    start = now()
    print("Collecting information from BAM file", file=sys.stderr)
    raw_barcodes = read_raw_barcodes(args.input)
    print(f"Found {len(raw_barcodes)} unique reads in {elapsed()}", file=sys.stderr)

    start = now()
    print(f"Analyzing barcodes", file=sys.stderr)
    barcode_mapping = collapse_barcodes(raw_barcodes, args)
    print(f"Corrected {len(barcode_mapping)} barcodes in {elapsed()}", file=sys.stderr)

    start = now()
    print(f"Creating output file", file=sys.stderr)
    convert_barcodes(args.input, args.output, barcode_mapping)
    print(f"Output file written in {elapsed()}", file=sys.stderr)

    if args.mappings:
        write_mappings(args.mappings, raw_barcodes, barcode_mapping)


def collapse_barcodes(raw_barcodes, args):
    """Run THE algorithm on all raw_barcodes

    Parameters
    ----------
    raw_barcodes: Set[str]
        iterable of raw_barcodes to collapse
    args: Parameters
        passed to THE algorithm

    Returns
    -------
    Dict[str, str]
        dict mapping of raw_barcode to collapsed barcode
    """
    barcodes = {}

    if args.perl_diff:
        # when args.perl_diff is set, all required get_len_differences
        # are computed at once, so only a single call to perl is required
        precompute_get_len_differences_using_perl(raw_barcodes, args)

    for i, raw_barcode in enumerate(raw_barcodes, 1):
        if i % 10000 == 0:
            print(f"\r{i}", end="", file=sys.stderr)
        try:
            barcodes[raw_barcode] = collapse_barcode(raw_barcode, args)
        except InvalidBarcode:
            pass
    print(f"\r{i}", file=sys.stderr)

    return barcodes


# THE algorithm


def collapse_barcode(raw_barcode, args):
    """Extract bc1 and bc2 from raw_barcode

    Parameters
    ----------
    raw_barcode: str
        barcode read from file
    args: Parameters
        namespace containing attributes barcodes1 and barcodes2
    Returns
    -------
    Optional[Tuple[str, str]]
        (raw_barcode, decoded_value) when found, None otherwise

    """
    bc_length = args.bc1_length
    if args.fixed_part1 in raw_barcode:
        # -- If it finds the first fixed part, looks at the first barcode next to it
        i = raw_barcode.index(args.fixed_part1)
        bc1_candidate = raw_barcode[:i][-bc_length:]
        bc1 = find_barcode(bc1_candidate, args.barcodes1, args.bc1_length, args)
    else:
        if args.both_fixed_parts_required:
            raise InvalidBarcode
        # fixed part not found, assume bc1 is in first bc_length characters
        bc1_candidate = raw_barcode[:bc_length]
        bc1 = find_barcode(bc1_candidate, args.barcodes1, args.bc1_length, args)
        # figure out how many characters were added/removed from fixed_part1
        fp1_candidate = raw_barcode[bc_length:bc_length+len(args.fixed_part1)]
        len_diff = get_len_difference(args.fixed_part1, fp1_candidate)
        i = bc_length + len_diff

    # discard part consumed by reading bc1
    tmp_barcode = raw_barcode[i + len(args.fixed_part1):][:args.bc2_length + len(args.fixed_part2) + 2]
    bc_length = args.bc2_length

    if args.fixed_part2 in tmp_barcode:
        # almost the same as for part1. one line is commented out
        # for comparability with perl script. If uncommented, this block
        # is the same as one above
        j = tmp_barcode.index(args.fixed_part2)
        bc2_candidate = tmp_barcode[:j]
        bc2 = find_barcode(bc2_candidate, args.barcodes2, args.bc2_length, args)
    else:
        if args.both_fixed_parts_required \
                or (args.one_fixed_part_required
                    and args.fixed_part1 not in raw_barcode):
            raise InvalidBarcode
        if len(tmp_barcode) < bc_length:
            tmp_barcode += "?" * (bc_length - len(tmp_barcode))
        bc2_candidate = tmp_barcode[:bc_length]
        bc2 = find_barcode(bc2_candidate, args.barcodes2, args.bc2_length, args)

    return bc1 + bc2


# dict for caching find_barcode results
cache = {}


def find_barcode(barcode, valid_barcodes, expected_length, args):
    """Find the barcode in the set of valid_barcodes

    If barcode has correct length, hamming distance is used
    Otherwise, levenshtein distance is used.

    If no candidates fit, or there are more than one equivalent candidates
    function raises InvalidBarcode

    Parameters
    ----------
    barcode: str
        barcode read from the file
    valid_barcodes: List[str]
        list of valid barcodes
    expected_length: int
        expected length of a valid barcode
    args: Parameters
        optional parameters
    Returns
    -------
    str
        closest matching candidate
    """

    # if barcode is already valid, return it immediately
    if barcode in valid_barcodes:
        return barcode

    # check if result has already been computed
    cache_key = (barcode, valid_barcodes)
    result = cache.get(cache_key)
    if result is None:
        # decide which distance to use based on barcode length
        if len(barcode) == expected_length:
            distance_func = hamming
            max_distance = args.max_distance
        else:
            distance_func = levenshtein
            # TODO: use abs to allow insertions
            max_distance = expected_length - len(barcode)

        # collect all candidates closer than max_distance
        candidates_at_distance = defaultdict(list)
        for candidate in valid_barcodes:
            d = distance_func(barcode, candidate)
            if d > max_distance:
                continue
            candidates_at_distance[d].append(candidate)

        # find the closest candidate
        for d in range(1, max_distance+1):
            if d in candidates_at_distance:
                if len(candidates_at_distance[d]) == 1:
                    # exactly one candidate, return it
                    result = candidates_at_distance[d][0]
                else:
                    # multiple candidates, do not select any
                    result = InvalidBarcode
                break
        else:
            # no barcode found at allowed distance
            result = InvalidBarcode

        # store result in cache
        cache[cache_key] = result

    if result is InvalidBarcode:
        raise result
    else:
        return result


get_len_difference_cache = {}


def get_len_difference(original, altered):
    """Return number of characters added/removed to/from `original`

    Assumes a "simple" transformation was performed.

    Parameters
    ----------
    original: str
        original sequence
    altered: str
        transformed sequence

    Returns
    -------
    int
        number of characters added to the original sequence (negative for deletions)
    """

    if altered in get_len_difference_cache:
        return get_len_difference_cache[altered]

    if args.perl_diff:
        # If args.perl_diff, all results should have been precomputed
        raise KeyError(f"Did not find results for {altered} in cache")

    if hamming(original, altered) < 2:
        # treat all replacements of 1 character as mutations that do not
        # change the length of the fixed part
        return 0

    # Use SequenceMatches to get a list of operations needed to transform
    # old into new
    t = difflib.SequenceMatcher(None, original, altered).get_opcodes()

    i = len([o for o, *_ in t if o == "insert"])
    d = len([o for o, *_ in t if o == "delete"])
    r = len([o for o, *_ in t if o == "replace"])
    if i == d == 1 and r == 0:
        # If there is exactly one insertion and deletion, find out which
        # occurs first and return how many characters it adds/removes
        for op, i1, i2, j1, j2 in t:
            if op == "delete":
                return i1-i2
            if op == "insert":
                return j2-j1
    if i == d == 0 and r == 1:
        # If there is only one replacement, return 0,
        # since the length of the sequence did not change
        return 0

    # If the change is more complex, return 0 (which is what perl script does).
    # A better approach would be to abort (raise InvalidBarcode) to avoid false positives
    return 0


def precompute_get_len_differences_using_perl(barcodes, args):
    """Compute get_len_differences using perl diff and store them
     in get_len_difference_cache

    Parameters
    ----------
    barcodes: List[str]
        iterable of all raw barcodes
    args: Parameters
    """
    # diffs will be needed for barcodes where fixed_part1 cannot be found
    fixed_part1_candidates = [
        barcode[args.bc1_length:][:len(args.fixed_part1)]
        for barcode in barcodes
        if args.fixed_part1 not in barcode
    ]

    # call perl to compute diffs
    diff = perl_diff(args.fixed_part1, fixed_part1_candidates)

    # compute differences and populate the cache
    for part1_candidate in fixed_part1_candidates:
        distance1, distance2 = diff[part1_candidate]
        count1 = distance1.count("[")
        pos1 = distance1.find("[")
        num_pos1 = distance1.find("]") - pos1 - 1

        count2 = distance2.count("{")
        pos2 = distance2.find("{")
        num_pos2 = distance2.find("}") - pos2 - 1

        d = 0
        if count1 == count2 == 1:
            if pos1 > pos2:
                d = num_pos1
            elif pos2 > pos1:
                d = -num_pos2

        get_len_difference_cache[part1_candidate] = d


def perl_diff(original, alternatives):
    """Compute diffs between original and each of the alternatives

    Parameters
    ----------
    original: str
    alternatives: List[str]

    Returns
    -------
    Dict[str, Tuple[str, str]]
        dict mapping alternative to result of perl diff
    """
    program = f"""
    use String::Diff 'diff';
    while(<>) {{
        chomp($_);
        ($distance_1, $distance_2) = diff("{original}", $_);
        print("$distance_1,$distance_2\n");
    }}
    """
    cmd = ["perl", "-e", program]
    results = run(cmd, input="\n".join(alternatives), stdout=PIPE, encoding="utf-8")

    diff = {}
    for input, result in zip(alternatives, results.stdout.split("\n")):
        result = result.strip().split(",")
        diff[input] = result
    return diff


class InvalidBarcode(Exception):
    """Exception raised when part of barcode cannot be found"""


# I/O functions

def read_bam_file(filename):
    """Open `filename` using samtools and read its output

    Parameters
    ----------
    filename: str
        name of the file to open

    Returns
    -------
    Generator[str]
        generator of lines from file
    """
    cmd = ["samtools", "view", "-h", filename]
    with Popen(cmd, stdout=PIPE, encoding="utf-8") as proc:
        for line in proc.stdout:
            yield line
        if proc.wait() != 0:
            raise CalledProcessError(proc.returncode, cmd)


def read_raw_barcodes(filename):
    """Read bam file, extracting RG:Z: part from non-header lines

    Parameters
    ----------
    filename: str
        path to the bam file

    Returns
    -------
    set[str]
    """
    raw_barcodes = set()
    for line in read_bam_file(filename):
        if line.startswith("@"):
            continue
        raw_barcodes.add(RG_RE.search(line).group(1))
    return raw_barcodes


def read_barcode_file(filename):
    """Read barcodes from file

    Parameters
    ----------
    filename: str
        path to the file containing valid bcs

    Returns
    -------
    frozenset[str]
    """
    with open(filename, "r") as f:
        return frozenset(bc.strip() for bc in f)


def convert_barcodes(input_bam, output_bam, barcode_mapping):
    """Copy input_bam to output_bam, converting barcodes using barcode_mapping

    Parameters
    ----------
    input_bam: str
        path to the input bam file
    output_bam: str
        path to the output bam file
    barcode_mapping: Dict[str, str]
        dict mapping raw barcodes to their collapsed versions
    """
    cmd = ["samtools", "view", "-b", "-@", "4", "-o", output_bam]
    with Popen(cmd, stdin=PIPE, encoding="utf-8") as proc:
        f = proc.stdin

        bam = read_bam_file(input_bam)
        header_written = False
        for line in bam:
            # Copy (non @RG) headers
            if line.startswith("@"):
                if not line.startswith("@RG"):
                    f.write(line)
                continue

            # Insert @RG lines before first "content" line
            if not header_written:
                for bc_RG_ID in sorted(set(barcode_mapping.values())):
                    f.write(f"@RG\tID:{bc_RG_ID}\tSM:{bc_RG_ID}\n")
                header_written = True

            # Copy "content" lines if barcode can be mapped
            match = RG_RE.search(line)
            raw_barcode = match.group(1)
            new_barcode = barcode_mapping.get(raw_barcode)
            if new_barcode is not None:
                f.write(line.replace(f"RG:Z:{raw_barcode}", f"RG:Z:{new_barcode}"))


def write_mappings(filename, raw_barcodes, barcode_mapping):
    """Write all raw_barcodes and their mappings into a file

    Parameters
    ----------
    filename: str
        path to the file where mappings should be written to
    raw_barcodes: frozenset(str)
        set containing all raw barcodes
    barcode_mapping: Dict[str, str]
        dict mapping raw barcodes to their collapsed versions
    """
    with open(filename, "w") as f:
        for raw_barcode in sorted(raw_barcodes):
            f.write(f"{raw_barcode}->{barcode_mapping.get(raw_barcode)}\n")


# Command line argument parsing


class Parameters(Namespace):
    input = None
    output = None
    mappings = None

    max_distance = 2

    barcodes1 = None
    barcodes2 = None
    bc1_length = None
    bc2_length = None
    fixed_part1 = "AGTACGTACGAGTC"
    fixed_part2 = "GTACTCGCAGTAGTC"

    one_fixed_part_required = False
    both_fixed_parts_required = False

    perl_diff = False

    @staticmethod
    def from_args():
        """Read parameters from command line args

        Returns
        -------
        Parameters
        """
        parser = ArgumentParser(description="Collapse Barcode Tags")
        parser.add_argument(
            "-i", "--input", required=True,
            help="input file name",
        )
        parser.add_argument(
            "-o", "--output", required=True,
            help="output file name",
        )
        parser.add_argument(
            "-d", "--max-distance", default=Parameters.max_distance, type=int,
            help="max hamming distance to consider two barcodes equal "
                 "[default: 2]",
        )

        chemistry_options = parser.add_argument_group("chemistry options")
        chemistry_options.add_argument(
            "-b", "--barcodes1", required=True, type=read_barcode_file,
            help="filename of file with valid barcodes for part1",
        )
        chemistry_options.add_argument(
            "-c", "--barcodes2", required=True, type=read_barcode_file,
            help="filename of file with valid barcodes for part2",
        )
        chemistry_options.add_argument(
            "--bc1-length",
            help="length of the bc1 segment [default: guessed from barcodes1]",
        )
        chemistry_options.add_argument(
            "--bc2-length",
            help="length of the bc2 segment [default: guessed from barcodes2]",
        )
        chemistry_options.add_argument(
            "-fp1", "--fixed-part1", default=Parameters.fixed_part1,
            help="fixed part after bc1 in sequence [default: AGTACGTACGAGTC]",
        )
        chemistry_options.add_argument(
            "-fp2", "--fixed-part2", default=Parameters.fixed_part2,
            help="fixed part after bc2 in sequence [default: GTACTCGCAGTAGTC]"
        )

        debug_options = parser.add_argument_group("debug options")
        debug_options.add_argument(
            "--one-fixed-part-required", action="store_true",
            help="ignore barcode if both of the fixed parts are missing/modified",
        )
        debug_options.add_argument(
            "--both-fixed-parts-required", action="store_true",
            help="ignore barcode if any of the fixed parts is missing/modified",
        )
        debug_options.add_argument(
            "--mappings",
            help="when specified, barcode mapping will be written to this file"
        )
        debug_options.add_argument(
            "--perl-diff", action="store_true",
            help="use perl's diff in get_len_difference, "
                 "use when same results as in perl script are required."
        )
        args = parser.parse_args(namespace=Parameters())

        if args.bc1_length is None:
            args.bc1_length = len(next(iter(args.barcodes1)))
        else:
            assert all(len(bc) == args.bc1_length for bc in args.barcodes1)
        if args.bc2_length is None:
            args.bc2_length = len(next(iter(args.barcodes2)))
        else:
            assert all(len(bc) == args.bc2_length for bc in args.barcodes2)
        return args


if __name__ == "__main__":
    args = Parameters.from_args()
    main(args)
