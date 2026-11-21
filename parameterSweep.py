"""
Parameter sweep over the visit settings on small random chains, meant for overnight runs.

Usage:
    python parameterSweep.py [--joints 4] [--seeds 8] [--first-seed 0] [--settings LABEL,...] [--no-checkpoints]
    python parameterSweep.py --resume "Trials Before Experiments/<sweep folder>"     # continue an interrupted sweep
    python parameterSweep.py --summarize "Trials Before Experiments/<sweep folder>"  # print the comparison table

Every combination of the settings in SETTINGS_GRID is run on each random chain (generated as in
randomChain.py, with collision-free configurations), one seed at a time, so stopping early still leaves
complete comparisons for the finished seeds. Each run uses the same traversal (DFS outward, one pass per
entry in tree.Joints) and resets np.random to the chain's seed first, so runs on the same chain differ only
in their settings.

Output, in "Trials Before Experiments/<timestamp>_Sweep_Joints<N>_Seeds<a>-<b>/":
    run_info.json            git commit, versions, and the sweep's settings
    summary.csv              one row per run: settings, start/final total length, sum of squares, time,
                             retries, reverts, final collisions, and any error
    errors.log               tracebacks of runs that raised
    seed<k>/start.tree, seed<k>/configs.json
    seed<k>/<settings label>/visits.csv, final.tree, checkpoints/, collision_penalty_log.txt
"""
import argparse
import copy
import csv
import itertools
import json
import os
import time
import traceback
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
from optimizationFunctions import *
from KinematicChain import *

SETTINGS_GRID = {
    'waypointSampling': ['box', 'multiscale'],
    'straightThroughGuess': [False, True],
    'alignmentWeight': [10, 0],
    'power': [2, 1],
}
TRAVERSAL = dict(traversal="dfs", direction="outward", orderBy="longest", repeatTraversal="n",
                 childFraction=1, guarantee=True)
CUBE_SIZE = 10
NUM_EXTRA_CONFIGS = 4  # collision-free random configurations in addition to neutral, as in randomChain.py

SUMMARY_FIELDS = ['seed', 'label'] + list(SETTINGS_GRID) + \
                 ['start_length', 'final_length', 'lower_bound', 'start_sum_squared', 'final_sum_squared',
                  'time_s', 'visits', 'retried', 'reverted', 'final_colliding_configs', 'error']

def settingsLabel(settings):
    return "_".join([settings['waypointSampling'],
                     "straight" if settings['straightThroughGuess'] else "nostraight",
                     f"align{settings['alignmentWeight']:g}",
                     f"power{settings['power']:g}"])

def allSettings():
    keys = list(SETTINGS_GRID)
    return [dict(zip(keys, values)) for values in itertools.product(*SETTINGS_GRID.values())]

def generateRandomChain(nJoints, r=1):
    # Same construction as randomChain.py
    bounds = (-CUBE_SIZE/2, CUBE_SIZE/2)
    poses = [SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds) for _ in range(nJoints)]
    jointConstructors = [
        lambda pose: TransverseRevolute(r=r, Pose=pose, minAngle=-np.pi/2, maxAngle=np.pi/2),
        lambda pose: CoaxialRevolute(r=r, Pose=pose, neutralLength=2, minAngle=-np.pi/2, maxAngle=np.pi/2),
        lambda pose: Prismatic(r=r, neutralLength=3, minLength=1, maxLength=5, Pose=pose)
    ]
    chain = KinematicChain(np.random.choice(jointConstructors)(poses[0]), gimbal=True)
    for i in range(1, nJoints):
        chain.appendGeneralizedGimbal(np.random.choice(jointConstructors)(poses[i]))
    return chain

def findCollisionFreeConfigs(tree, maxAttempts=100, numNeeded=NUM_EXTRA_CONFIGS):
    # Same as randomChain.py's find_collision_free_configs: neutral plus up to numNeeded random ones
    configs = [[0] * len(tree.Joints)]
    for _ in range(maxAttempts):
        if len(configs) >= numNeeded + 1:
            break
        config = tree.randomConfiguration(realJointsOnly=False).tolist()
        test = copy.deepcopy(tree)
        for i in range(len(test.Joints)):
            if not isWaypoint(test.Joints[i]):
                test.setJointState(i, config[i])
            test.Joints[i].recomputeCollisionCapsules()
        if test.detectCollisions(debug=False) == 0:
            configs.append(config)
    return configs

def prepareSeed(sweepDir, seed, nJoints):
    seedDir = os.path.join(sweepDir, f"seed{seed}")
    treePath = os.path.join(seedDir, "start.tree")
    configsPath = os.path.join(seedDir, "configs.json")
    if os.path.exists(treePath) and os.path.exists(configsPath):
        return seedDir, loadKinematicTree(treePath), json.load(open(configsPath))
    os.makedirs(seedDir, exist_ok=True)
    np.random.seed(seed)
    chain = generateRandomChain(nJoints)
    configs = findCollisionFreeConfigs(chain)
    with open(treePath, 'w') as f:
        f.write(repr(chain))
    with open(configsPath, 'w') as f:
        json.dump(configs, f, indent=2)
    # use the reloaded tree, so runs and resumed runs start from exactly the same object state
    return seedDir, loadKinematicTree(treePath), configs

def appendSummary(sweepDir, row):
    path = os.path.join(sweepDir, 'summary.csv')
    isNew = not os.path.exists(path)
    with open(path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, restval='', extrasaction='ignore')
        if isNew:
            writer.writeheader()
        writer.writerow(row)

def finishedRuns(sweepDir):
    path = os.path.join(sweepDir, 'summary.csv')
    if not os.path.exists(path):
        return set()
    with open(path) as f:
        return {(int(r['seed']), r['label']) for r in csv.DictReader(f)}

def runOne(sweepDir, seedDir, seed, chain, configs, settings, saveCheckpoints):
    label = settingsLabel(settings)
    runDir = os.path.join(seedDir, label)
    os.makedirs(runDir, exist_ok=True)
    setCollisionPenaltyLogFile(os.path.join(runDir, 'collision_penalty_log.txt'))
    row = dict(seed=seed, label=label, **settings, start_length=chain.totalLength(),
               lower_bound=chain.totalLengthLowerBound(), start_sum_squared=sumSquaredLinkLength(chain))
    start = time.time()
    try:
        np.random.seed(seed)
        tree, times, lengths = optimizeTree(copy.deepcopy(chain), configurations=configs, showSteps=False,
                                            parallelize=True, evaluate=True, verbose=False, directory=runDir,
                                            save_checkpoints=saveCheckpoints, **TRAVERSAL, **settings)
        with open(os.path.join(runDir, 'visits.csv')) as f:
            visits = [r for r in csv.DictReader(f) if r['visit'] != '0']
        row.update(final_length=tree.totalLength(), final_sum_squared=sumSquaredLinkLength(tree),
                   visits=len(visits), retried=sum(r['retried'] == 'True' for r in visits),
                   reverted=sum(r['reverted'] == 'True' for r in visits),
                   final_colliding_configs=sum(len(collidingConfigurations(tree, i, configs))
                                               for i in range(len(tree.Joints))))
    except Exception as e:
        row['error'] = f"{type(e).__name__}: {e}"
        with open(os.path.join(sweepDir, 'errors.log'), 'a') as f:
            f.write(f"seed {seed}, {label}:\n{traceback.format_exc()}\n")
    row['time_s'] = time.time() - start
    appendSummary(sweepDir, row)
    return row

def summarize(sweepDir):
    with open(os.path.join(sweepDir, 'summary.csv')) as f:
        rows = [r for r in csv.DictReader(f)]
    ok = [r for r in rows if not r['error']]
    failed = [r for r in rows if r['error']]
    seeds = sorted({int(r['seed']) for r in ok})
    baseline = settingsLabel(allSettings()[0])
    base = {int(r['seed']): float(r['final_length']) for r in ok if r['label'] == baseline}
    print(f"{sweepDir}\n{len(ok)} finished runs on seeds {seeds}, {len(failed)} failed\n")
    print(f"Final total length relative to the baseline ({baseline}) on the same chain; lower is better.")
    print(f"{'settings':45} {'runs':>4} {'mean ratio':>10} {'best':>6} {'worst':>6} {'wins':>5} "
          f"{'mean final':>10} {'retried':>7} {'reverted':>8} {'mean time':>9}")
    results = []
    for label in [settingsLabel(s) for s in allSettings()]:
        runs = [r for r in ok if r['label'] == label and int(r['seed']) in base]
        if not runs:
            continue
        ratios = [float(r['final_length']) / base[int(r['seed'])] for r in runs]
        results.append((np.mean(ratios), label, runs, ratios))
    for meanRatio, label, runs, ratios in sorted(results):
        wins = sum(ratio < 1 - 1e-9 for ratio in ratios)
        print(f"{label:45} {len(runs):>4} {meanRatio:>10.3f} {min(ratios):>6.3f} {max(ratios):>6.3f} {wins:>5} "
              f"{np.mean([float(r['final_length']) for r in runs]):>10.1f} "
              f"{sum(int(r['retried']) for r in runs):>7} {sum(int(r['reverted']) for r in runs):>8} "
              f"{np.mean([float(r['time_s']) for r in runs]) / 60:>8.1f}m")
    if failed:
        print("\nFailed runs (see errors.log):")
        for r in failed:
            print(f"  seed {r['seed']} {r['label']}: {r['error']}")

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--joints', type=int, default=4, help="joints per random chain (default 4)")
    parser.add_argument('--seeds', type=int, default=8, help="number of random chains (default 8)")
    parser.add_argument('--first-seed', type=int, default=0, help="seed of the first chain (default 0)")
    parser.add_argument('--settings', help="comma-separated settings labels to run (default: the whole grid)")
    parser.add_argument('--no-checkpoints', action='store_true', help="don't save a checkpoint after every visit")
    parser.add_argument('--resume', metavar='SWEEP_DIR', help="continue an interrupted sweep in this folder")
    parser.add_argument('--summarize', metavar='SWEEP_DIR', help="print the comparison table for a sweep and exit")
    args = parser.parse_args()

    if args.summarize:
        summarize(args.summarize)
        return

    if args.resume:
        sweepDir = args.resume
        with open(os.path.join(sweepDir, 'run_info.json')) as f:
            sweep = json.load(f)['settings']
    else:
        sweep = dict(script='parameterSweep.py', joints=args.joints,
                     seeds=list(range(args.first_seed, args.first_seed + args.seeds)),
                     settings=args.settings.split(',') if args.settings else [settingsLabel(s) for s in allSettings()],
                     grid=SETTINGS_GRID, traversal=TRAVERSAL, cubeSize=CUBE_SIZE,
                     numExtraConfigs=NUM_EXTRA_CONFIGS, saveCheckpoints=not args.no_checkpoints)
        timestamp = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
        sweepDir = os.path.join("Trials Before Experiments",
                                f"{timestamp}_Sweep_Joints{args.joints}_Seeds{sweep['seeds'][0]}-{sweep['seeds'][-1]}")
        os.makedirs(sweepDir, exist_ok=True)
        writeRunInfo(sweepDir, sweep)

    settingsList = [s for s in allSettings() if settingsLabel(s) in sweep['settings']]
    unknown = set(sweep['settings']) - {settingsLabel(s) for s in allSettings()}
    if unknown:
        raise SystemExit(f"Unknown settings labels {sorted(unknown)}; valid: {[settingsLabel(s) for s in allSettings()]}")
    done = finishedRuns(sweepDir)
    total = len(sweep['seeds']) * len(settingsList)
    print(f"Sweep folder: {sweepDir}\n{total} runs ({len(sweep['seeds'])} chains x {len(settingsList)} settings), "
          f"{len(done)} already finished\n", flush=True)

    runTimes = []
    count = len(done)
    for seed in sweep['seeds']:
        seedDir, chain, configs = prepareSeed(sweepDir, seed, sweep['joints'])
        for settings in settingsList:
            label = settingsLabel(settings)
            if (seed, label) in done:
                continue
            row = runOne(sweepDir, seedDir, seed, chain, configs, settings, sweep['saveCheckpoints'])
            count += 1
            runTimes.append(row['time_s'])
            remaining = (total - count) * np.mean(runTimes)
            result = f"error: {row['error']}" if row.get('error') else \
                     f"{row['start_length']:.1f} -> {row['final_length']:.1f}, retried {row['retried']}, reverted {row['reverted']}"
            print(f"[{datetime.now():%H:%M:%S}] {count}/{total} seed {seed} {label}: {result} "
                  f"({row['time_s'] / 60:.1f} min; about {remaining / 3600:.1f} h left)", flush=True)
    print()
    summarize(sweepDir)

if __name__ == "__main__":
    main()
