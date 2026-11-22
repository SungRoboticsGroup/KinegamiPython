"""Usage: python debugExamples/bunching/surveyBunching.py <tree file> <configs json> <output json>

For each waypoint of a tree, search its placements much harder than the optimizer does, to find how much
shorter its adjacent links could be. Reports both the optimizer's waypoint objective (linkLoss + 10 *
misalignment with the parent) and plain linkLoss, and saves the best moves found for each to the output json.

This produced the reference placements in examples.json (run on the final tree of the 8-joint seed-44 chain,
2026-11-21; about 90 seconds per waypoint with 16 processes)."""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from optimizationFunctions import *

treePath, configsPath, outPath = sys.argv[1], sys.argv[2], sys.argv[3]
N_RANDOM = 3000
tree = loadKinematicTree(treePath)
configs = json.load(open(configsPath))

def aligned(w):
    wp, parent = tree.Joints[w], tree.Joints[tree.Parents[w]]
    alignTransform = parent.DistalDubinsFrame() * wp.ProximalDubinsFrame().inv()
    return alignTransform.R @ wp.Pose.R

def evaluate(args):
    w, A = args
    c = tree.copyAbbreviatedSelf()
    if not c.transformJoint(w, SE3(A, check=False), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
        return (np.inf, np.inf)
    L = linkLoss(c, w)
    mis = rotationAngle(c.Joints[w].Pose.R.T @ aligned(w))
    return (L, mis)

def sampleMoves(w, rng):
    d = localLinkScale(tree, w)
    moves = []
    for k in range(N_RANDOM):
        R = Rotation.random(random_state=rng).as_matrix()
        if k < N_RANDOM // 3:
            t = np.zeros(3)                       # rotate in place
        elif k < 2 * N_RANDOM // 3:
            t = rng.uniform(-0.5 * d, 0.5 * d, 3)  # small moves
        else:
            t = rng.uniform(-2 * d, 2 * d, 3)      # same range as the optimizer's bounds
        A = np.eye(4); A[:3, :3] = R; A[:3, 3] = t
        moves.append(A)
    return moves

def scoreObjective(L, mis):
    return L + 10 * mis

def scoreLinkLoss(L, mis):
    return L

SCORES = {'objective': scoreObjective, 'linkLoss': scoreLinkLoss}

def refine(w, A0, scoreName):
    # Nelder-Mead on a small local perturbation of the move A0
    score = SCORES[scoreName]
    def f(p):
        A = A0 @ SE3.Trans(p[:3]).A @ SE3.RPY(p[3:]).A
        L, mis = evaluate((w, A))
        return score(L, mis)
    res = minimize(f, np.zeros(6), method='Nelder-Mead', options={'maxiter': 400, 'xatol': 1e-3, 'fatol': 1e-4})
    return A0 @ SE3.Trans(res.x[:3]).A @ SE3.RPY(res.x[3:]).A, res.fun

if __name__ == '__main__':
    rng = np.random.default_rng(0)
    results = []
    with Pool(16) as pool:
        for w in [i for i in range(len(tree.Joints)) if isWaypoint(tree.Joints[i])]:
            start = time.time()
            L0, mis0 = evaluate((w, np.eye(4)))
            moves = sampleMoves(w, rng)
            evals = pool.map(evaluate, [(w, A) for A in moves])
            entry = {'index': w, 'incoming': tree.Links[w].path.length, 'scale': localLinkScale(tree, w),
                     'current': {'linkLoss': L0, 'misalignment': mis0, 'objective': L0 + 10 * mis0},
                     'infeasible': int(sum(np.isinf(e[0]) for e in evals))}
            for name, score in SCORES.items():
                order = np.argsort([score(*e) for e in evals])[:4]
                refined = pool.starmap(refine, [(w, moves[i], name) for i in order])
                A, best = min(refined, key=lambda r: r[1])
                L, mis = evaluate((w, A))
                moved = tree.copyAbbreviatedSelf()
                moved.transformJoint(w, SE3(A, check=False), propogate=False, safe=True, relative=True, recomputeBoundingBall=False)
                entry[name] = {'linkLoss': L, 'misalignment': mis, 'objective': L + 10 * mis,
                               'colliding': len(collidingConfigurations(moved, w, configs)),
                               'translation': float(np.linalg.norm(A[:3, 3])),
                               'rotation': float(rotationAngle(A[:3, :3])), 'move': A.tolist()}
            results.append(entry)
            c, o, l = entry['current'], entry['objective'], entry['linkLoss']
            print(f"W{w:2d} ({time.time()-start:4.0f}s): current linkLoss {c['linkLoss']:8.2f} mis {c['misalignment']:.2f} obj {c['objective']:8.2f} | "
                  f"best obj {o['objective']:8.2f} (linkLoss {o['linkLoss']:8.2f}, coll {o['colliding']}) | "
                  f"best linkLoss {l['linkLoss']:8.2f} (mis {l['misalignment']:.2f}, obj {l['objective']:8.2f}, coll {l['colliding']}, "
                  f"move {l['translation']:.2f}, rot {np.degrees(l['rotation']):.0f}deg)", flush=True)
    json.dump(results, open(outPath, 'w'), indent=1)
