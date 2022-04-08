import argparse
import functools
from operator import attrgetter
from typing import List
import heapq
import contextlib


import dnaio
from xopen import xopen

READ_OPENER = functools.partial(xopen, threads=1)


def fastq_to_sorted_chunks(fastq: str, output_prefix: str,
                           max_ram_records: int = 5_000_000) -> List[str]:
    sorted_chunks = []
    with dnaio.open(fastq, mode="r", opener=READ_OPENER) as fastq_reader:
        while True:
            filename = output_prefix + str(len(sorted_chunks))
            records = [record for i, record in
                       zip(range(max_ram_records), fastq_reader)]
            if len(records) == 0:
                break
            records.sort(key=attrgetter('sequence'))
            with xopen(filename, mode="wb", format="gz", threads=1, compresslevel=1
                       ) as sorted_chunk:
                for record in records:
                    sorted_chunk.write(record.fastq_bytes())
    return sorted_chunks


def sorted_chunks_to_sorted_fastq(sorted_chunks: List[str], output_file: str):
    with contextlib.ExitStack() as stack:
        files = [stack.enter_context(
                 dnaio.open(sorted_chunk, mode="r", opener=READ_OPENER))
                 for sorted_chunk in sorted_chunks]
        with heapq.merge(files, key=attrgetter('sequence')) as merged:
            with xopen(output_file, mode="wb", compresslevel=1, threads=1
                       ) as output:
                for record in merged:
                    output.write(record.fastq_bytes())


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", metavar="INPUT", help="Input FASTQ file")
    parser.add_argument("-o", "--output", help="Output FASTQ file")
    return parser


def main():
    args = argument_parser().parse_args()
    sorted_chunks = fastq_to_sorted_chunks(args.input, args.output)
    sorted_chunks_to_sorted_fastq(sorted_chunks, args.output)


if __name__ == "__main__":
    main()

