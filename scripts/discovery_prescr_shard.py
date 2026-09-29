"""Deterministic MHCflurry prescreen shard, matching discovery_funnel.py chunk boundaries.

A run selects half-open chunk indices [start, end). Existing chunks are never
recomputed. Input order comes from external sort -u of frozen funnel_self.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys

from discovery_funnel import PANEL, CHUNK, D, _mhcflurry

INPUT_SHA = 'cf9b52a3a997e1c07a76fa5168375aad4638a5317e6d054dfa6deaf76bf82064'
N_PEPTIDES = 8212079


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', type=int, required=True)
    ap.add_argument('--end', type=int, required=True)
    args = ap.parse_args()
    per_chunk = CHUNK // len(PANEL)
    total = (N_PEPTIDES + per_chunk - 1) // per_chunk
    assert 0 <= args.start < args.end <= total, (args.start, args.end, total)
    input_file = f'{D}/_uniq_self.txt.gz'
    if not os.path.exists(input_file):
        # Same sorted unique peptide operation as phase_prescr, using the
        # already frozen self-filtered records and POSIX sort under C locale.
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            unsorted = os.path.join(td, 'peptides.txt')
            with gzip.open(f'{D}/funnel_self.jsonl.gz', 'rt') as fin, open(unsorted, 'w') as out:
                for line in fin:
                    for pep in json.loads(line)['peptides']:
                        out.write(pep + '\n')
            with open(input_file, 'wb') as out:
                sorter = subprocess.Popen(['sort', '-u', unsorted], env={**os.environ, 'LC_ALL': 'C'}, stdout=subprocess.PIPE)
                zipper = subprocess.Popen(['gzip', '-n'], stdin=sorter.stdout, stdout=out)
                sorter.stdout.close()
                assert zipper.wait() == 0 and sorter.wait() == 0
    sha = hashlib.sha256()
    count = 0
    peptides = []
    with gzip.open(input_file, 'rb') as fin:
        for raw in fin:
            sha.update(raw)
            if args.start * per_chunk <= count < args.end * per_chunk:
                peptides.append(raw.decode().strip())
            count += 1
    assert count == N_PEPTIDES and sha.hexdigest() == INPUT_SHA, (count, sha.hexdigest())
    from mhcflurry import Class1AffinityPredictor
    pred = Class1AffinityPredictor.load()
    os.makedirs(f'{D}/mhc_chunks', exist_ok=True)
    for chunk in range(args.start, args.end):
        path = f'{D}/mhc_chunks/chunk_{chunk:05d}.json.gz'
        if os.path.exists(path):
            continue
        lo = (chunk - args.start) * per_chunk
        group = peptides[lo:lo + per_chunk]
        assert len(group) == min(per_chunk, N_PEPTIDES - chunk * per_chunk)
        pairs = [(p, a) for p in group for a in PANEL]
        output = _mhcflurry(pred, pairs)
        assert len(output) == len(pairs)
        assert all(tuple(row[:2]) == pair for row, pair in zip(output, pairs))
        tmp = path + '.partial'
        with gzip.open(tmp, 'wt') as fh:
            json.dump(output, fh)
        os.replace(tmp, path)
        print(f'chunk {chunk}/{total - 1} {len(output)} pairs', flush=True)
        if (chunk + 1) % 20 == 0:
            os.execv(sys.executable, [sys.executable] + sys.argv)


if __name__ == '__main__':
    main()
