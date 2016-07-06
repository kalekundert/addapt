#!/usr/bin/env python3

import numpy as np
import pandas as pd

def load_trajectory(tsv_path):
    return pd.read_table(tsv_path or 'logs/mh.tsv')

def load_trajectories(tsv_paths):
    return [load_trajectory(x) for x in (tsv_paths or [None])]

def pick_best_seqs(tsv_paths, window_size):
    best_seqs = []

    for traj in load_trajectories(tsv_paths):
        best_scores = pick_best_scores(traj['current_score'], window_size)
        best_seqs += [extract_pick(traj, i) for i in best_scores]

    return best_seqs

def pick_best_scores(scores, window_size):
    peaks = []
    scores = np.array(scores)
    BLOCKED = -1

    while scores.max() != BLOCKED:
        i = scores.argmax()
        window = slice(
                max(i - window_size//2, 0),
                min(i + window_size//2, len(scores)))

        if not any(scores[window] == BLOCKED):
            peaks.append(i)

        scores[window] = BLOCKED

    return sorted(peaks)

def extract_pick(traj, i):
    class Pick: pass    # (no fold)
    pick = Pick()

    pick.i = i
    pick.traj = traj
    pick.score = traj['current_score'][i]

    pick.seq = ''
    for col in traj.columns:
        if col.startswith('current_domain'):
            domain_seq = traj[col][i]
            if isinstance(domain_seq, str):
                pick.seq += domain_seq

    return pick


def run_rnafold(seq, theo=False):
    from subprocess import Popen, PIPE

    # Construct the RNAfold command line.  Include the theophylline aptamer 
    # motif if requested.
    cmd = 'RNAfold', '--noPS', '--MEA'
    if theo:
        cmd += '--motif', '{},{},{}'.format(
            'GAUACCAGCCGAAAGGCCCUUGGCAGC',
            '(...((.(((....)))....))...)',
            -9.22) # kcal/mol

    # Run the command and return the resulting lines as a string.
    process = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(seq.encode())
    return stdout.decode().strip()

