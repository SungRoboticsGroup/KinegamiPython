/**
 * dubinsCSC.js — client-side 3-D Dubins CSC path solver + link geometry
 *
 * Ported from PathCSC.py, geometryHelpers.py, LinkCSC.py.
 *
 * SE3 is represented as { t:[x,y,z], x:[…], y:[…], z:[…] }
 * where x/y/z are the *columns* of the rotation matrix.
 */

// ── Vec3 math ─────────────────────────────────────────────────────────────────

const vadd   = (a, b) => [a[0]+b[0], a[1]+b[1], a[2]+b[2]];
const vsub   = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const vscale = (a, s) => [a[0]*s,    a[1]*s,    a[2]*s];
const vdot   = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
const vcross = (a, b) => [
  a[1]*b[2] - a[2]*b[1],
  a[2]*b[0] - a[0]*b[2],
  a[0]*b[1] - a[1]*b[0],
];
const vnorm  = a => Math.sqrt(vdot(a, a));
const vunit  = a => { const n = vnorm(a); return n > 0 ? vscale(a, 1/n) : [1, 0, 0]; };

/** Any unit vector perpendicular to u (matching null_space([u])[:,0] semantics). */
function anyPerp(u) {
  const ax = Math.abs(u[0]), ay = Math.abs(u[1]), az = Math.abs(u[2]);
  if (ax <= ay && ax <= az) return vunit(vcross(u, [1, 0, 0]));
  if (ay <= az)             return vunit(vcross(u, [0, 1, 0]));
  return                           vunit(vcross(u, [0, 0, 1]));
}

/** Unit vector orthogonal to both u and v (right-hand-rule; fallback when parallel). */
function unitNormalToBoth(u, v) {
  const c = vcross(u, v);
  const n = vnorm(c);
  if (n < 1e-8) return anyPerp(u);
  return vscale(c, 1 / n);
}

/** Signed angle (radians) from a to b measured around n. */
function signedAngle(a, b, n) {
  const au = vunit(a), bu = vunit(b);
  const nu = vnorm(n) > 0 ? vunit(n) : n;
  return Math.atan2(vdot(vcross(au, bu), nu), vdot(au, bu));
}

/** Wrap angle to [0, 2π). */
function wrapAngle(angle) {
  const TAU = 2 * Math.PI;
  return ((angle % TAU) + TAU) % TAU;
}

/** IEEE 754-style remainder — wraps to (−π, π]. Matches Python math.remainder. */
function mathRemainder(x, y) {
  return x - Math.round(x / y) * y;
}

// ── SE3 math ──────────────────────────────────────────────────────────────────
// Pose: { t:[x,y,z], x:[…], y:[…], z:[…] }
//   t   = translation
//   x/y/z = columns 0/1/2 of the rotation matrix

const se3id = () => ({ t: [0,0,0], x: [1,0,0], y: [0,1,0], z: [0,0,1] });
const se3Tx = dw  => ({ t: [dw,0,0], x: [1,0,0], y: [0,1,0], z: [0,0,1] });

/** Apply rotation part of pose to vector v:  R @ v = x*v[0] + y*v[1] + z*v[2]. */
function mat3apply(pose, v) {
  return [
    pose.x[0]*v[0] + pose.y[0]*v[1] + pose.z[0]*v[2],
    pose.x[1]*v[0] + pose.y[1]*v[1] + pose.z[1]*v[2],
    pose.x[2]*v[0] + pose.y[2]*v[1] + pose.z[2]*v[2],
  ];
}

/** Apply SE3 to a point:  R @ p + t. */
const se3apply = (pose, p) => vadd(mat3apply(pose, p), pose.t);

/** SE3 composition:  A @ B. */
function se3mul(A, B) {
  return {
    t: vadd(mat3apply(A, B.t), A.t),
    x: mat3apply(A, B.x),
    y: mat3apply(A, B.y),
    z: mat3apply(A, B.z),
  };
}

/** Pure-rotation SE3 from angle-axis (Rodrigues formula). */
function se3AngleAxis(angle, axis) {
  const c = Math.cos(angle), s = Math.sin(angle), omc = 1 - c;
  const [ux, uy, uz] = vunit(axis);
  return {
    t: [0, 0, 0],
    x: [c + omc*ux*ux,     omc*ux*uy + uz*s,  omc*ux*uz - uy*s],
    y: [omc*ux*uy - uz*s,  c + omc*uy*uy,     omc*uy*uz + ux*s],
    z: [omc*ux*uz + uy*s,  omc*uy*uz - ux*s,  c + omc*uz*uz   ],
  };
}

// ── Newton-Raphson fsolve (4×4) ───────────────────────────────────────────────

/** Gaussian elimination with partial pivoting.  Solves A x = b, returns x or null. */
function gaussElim(A, b) {
  const n = 4;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let r = col + 1; r < n; r++)
      if (Math.abs(M[r][col]) > Math.abs(M[maxRow][col])) maxRow = r;
    [M[col], M[maxRow]] = [M[maxRow], M[col]];
    if (Math.abs(M[col][col]) < 1e-14) return null;
    for (let r = col + 1; r < n; r++) {
      const f = M[r][col] / M[col][col];
      for (let j = col; j <= n; j++) M[r][j] -= f * M[col][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

/** Newton-Raphson solver for F(x)=0, x∈R^4. Matches scipy.optimize.fsolve behaviour. */
function fsolve4(F, x0, maxIter = 60, tol = 1e-9) {
  const h = 1e-7;
  let x = [...x0];
  for (let iter = 0; iter < maxIter; iter++) {
    const fx = F(x);
    if (vnorm(fx) < tol) break;
    // Numerical Jacobian — J[j] is column j = dF/dx_j
    const J = x.map((_, j) => {
      const xp = [...x]; xp[j] += h;
      return F(xp).map((v, i) => (v - fx[i]) / h);
    });
    // Build row-major matrix A where A[i][j] = J[j][i]
    const A = Array.from({ length: 4 }, (_, i) => J.map(col => col[i]));
    const dx = gaussElim(A, fx.map(v => -v));
    if (!dx) break;
    for (let j = 0; j < 4; j++) x[j] += dx[j];
  }
  return x;
}

// ── PathCSC construction ──────────────────────────────────────────────────────
//
// Ported from PathCSC.py.
// Solves for the straight segment t of a Circle-Straight-Circle Dubins path
// given start/end Dubins frames and circle-direction signs.

function makePathCSC(tDirMag, r, sp, sd, ep, ed, c1s, c2s) {
  const tDir  = tDirMag.slice(0, 3);
  const tUnit = vunit(tDir);
  const tMag  = Math.abs(tDirMag[3]);
  const t     = vscale(tUnit, tMag);

  // Circle 1 (start)
  const cn1       = unitNormalToBoth(tUnit, sd);
  const w1        = vscale(vcross(sd,    cn1), c1s);
  const y1        = vscale(vcross(tUnit, cn1), c1s);
  const cc1       = vadd(sp, vscale(w1, r));
  const turn1end  = vsub(cc1, vscale(y1, r));

  // Circle 2 (end)
  const cn2        = unitNormalToBoth(tUnit, ed);
  const w2         = vscale(vcross(ed,    cn2), c2s);
  const y2         = vscale(vcross(tUnit, cn2), c2s);
  const cc2        = vadd(ep, vscale(w2, r));
  const turn2start = vsub(cc2, vscale(y2, r));

  // Error: t should equal (turn2start − turn1end); also enforce |tDir|=1
  const tError = vsub(t, vsub(turn2start, turn1end));
  const error  = [...tError, vnorm(tDir) - 1];

  // Arc angles (wrapped to [0, 2π))
  const theta1 = wrapAngle(signedAngle(sd,    tUnit, vcross(sd,    w1)));
  const theta2 = wrapAngle(signedAngle(tUnit, ed,    vcross(ed,    w2)));
  const length = r * theta1 + tMag + r * theta2;

  return { tUnit, tMag, t, w1, w2, y1, y2, cc1, cc2,
           turn1end, turn2start, theta1, theta2, length,
           error, c1s, c2s };
}

// ── CSC solver ────────────────────────────────────────────────────────────────

function _solveOne(r, sp, sd, ep, ed, c1s, c2s, x0) {
  const sol = fsolve4(x => makePathCSC(x, r, sp, sd, ep, ed, c1s, c2s).error, x0);
  return makePathCSC(sol, r, sp, sd, ep, ed, c1s, c2s);
}

/**
 * Find the shortest valid CSC Dubins path.
 * prevSolution = { tUnit, tMag, c1s, c2s } for warm-starting (or null).
 * Returns a path object with all PathCSC fields.
 */
export function shortestCSC(r, sp, sd, ep, ed, prevSolution = null) {
  const EPS = 1e-6;

  // Degenerate case: start = end
  if (vnorm(vsub(sp, ep)) < EPS && vnorm(vsub(sd, ed)) < EPS)
    return makePathCSC([...sd, 0], r, sp, sd, ep, ed, 1, 1);

  // Default initial guess: straight line from sp to ep
  let t0 = vsub(ep, sp);
  if (vnorm(t0) < EPS) t0 = vscale(sd, r);
  const x0 = [...vunit(t0), vnorm(t0)];

  const isValid = p =>
    vnorm(p.error) <= 0.005 * r &&
    p.theta1 < Math.PI &&
    p.theta2 < Math.PI;

  // Warm-start: try previous sign combination first with previous solution as x0
  if (prevSolution) {
    const { tUnit: pu, tMag: pm, c1s: pc1, c2s: pc2 } = prevSolution;
    const p = _solveOne(r, sp, sd, ep, ed, pc1, pc2, [...pu, pm]);
    if (isValid(p)) return p;
  }

  // Full search: try all 4 sign combinations
  const combos  = [[1,1],[1,-1],[-1,1],[-1,-1]];
  const paths   = combos.map(([c1, c2]) => _solveOne(r, sp, sd, ep, ed, c1, c2, x0));
  const lengths = paths.map(p => isValid(p) ? p.length : Infinity);
  return paths[lengths.indexOf(Math.min(...lengths))];
}

// ── Elbow geometry ────────────────────────────────────────────────────────────
//
// Ported from geometryHelpers.Elbow.circleEllipseCircleQT
// Generates a toroidal section: StartCircle → MidEllipse → EndCircle

function elbowGeo(r, StartFrame, bendingAngle, rotAxisAngle, numSides) {
  const EPS = 0.0001;

  // Normalise — match Python Elbow.__init__ exactly
  rotAxisAngle = ((rotAxisAngle % (2*Math.PI)) + 2*Math.PI) % (2*Math.PI);
  bendingAngle = mathRemainder(bendingAngle, 2*Math.PI);
  if (bendingAngle < 0) {
    bendingAngle = -bendingAngle;
    rotAxisAngle = ((rotAxisAngle + Math.PI) % (2*Math.PI));
  }
  if (bendingAngle < EPS) return null;

  // rotAxisDirLocal = SO3.Rx(rotAxisAngle) @ [0,1,0] = [0, cos, sin]
  const rotAxisDirLocal = [0, Math.cos(rotAxisAngle), Math.sin(rotAxisAngle)];

  const dw              = r * Math.tan(bendingAngle / 2);
  const Forward         = se3Tx(dw);
  const Rotate          = se3AngleAxis(bendingAngle, rotAxisDirLocal);
  const Transformation  = se3mul(Forward, se3mul(Rotate, Forward));
  const EndFrame        = se3mul(StartFrame, Transformation);

  // Mid-plane: normal = x-col of (StartFrame @ HalfRotate)
  const HalfRotate     = se3AngleAxis(bendingAngle / 2, rotAxisDirLocal);
  const midPoint       = se3apply(StartFrame, [dw, 0, 0]);       // (StartFrame@Forward).t
  const midPlaneNormal = se3mul(StartFrame, HalfRotate).x;       // x-col

  // dot(-ahat, midPlaneNormal)  — denominator for plane intersection
  const alignment = vdot(vscale(StartFrame.x, -1), midPlaneNormal);

  // Build 3 rings of vertices: StartCircle (ring 0), MidEllipse (ring 1), EndCircle (ring 2)
  const n = numSides;
  const vertices = [];

  for (let ring = 0; ring < 3; ring++) {
    for (let j = 0; j < n; j++) {
      const angle = 2 * Math.PI * j / n;
      const u = r * Math.cos(angle);
      const v = r * Math.sin(angle);

      if (ring === 0) {
        // StartCircle:  StartFrame.t + u*bhat + v*chat
        vertices.push(vadd(StartFrame.t,
          vadd(vscale(StartFrame.y, u), vscale(StartFrame.z, v))));
      } else if (ring === 1) {
        // MidEllipse: project StartCircle along ahat onto mid-plane
        const sc  = vadd(StartFrame.t,
          vadd(vscale(StartFrame.y, u), vscale(StartFrame.z, v)));
        const ts  = vdot(vsub(sc, midPoint), midPlaneNormal) / alignment;
        vertices.push(vadd(sc, vscale(StartFrame.x, ts)));
      } else {
        // EndCircle:  EndFrame.t + u*endBhat + v*endChat
        vertices.push(vadd(EndFrame.t,
          vadd(vscale(EndFrame.y, u), vscale(EndFrame.z, v))));
      }
    }
  }

  // Faces: 2 triangles per quad between consecutive rings
  const faces = [];
  for (let ring = 0; ring < 2; ring++) {
    for (let j = 0; j < n; j++) {
      const nj = (j + 1) % n;
      const a = ring*n + j,   b = ring*n + nj;
      const c = (ring+1)*n + j, d = (ring+1)*n + nj;
      faces.push([a, b, c]);
      faces.push([c, b, d]);
    }
  }

  return { vertices, faces, endFrame: EndFrame };
}

// ── Compound elbow geometry ───────────────────────────────────────────────────
//
// Ported from geometryHelpers.CompoundElbow.circleEllipseCircleQT

function compoundElbowGeo(r, StartFrame, bendingAngle, rotAxisAngle, maxAngle, numSides) {
  const EPS = 0.0001;

  // Normalise — match Python CompoundElbow.__init__
  rotAxisAngle = ((rotAxisAngle % (2*Math.PI)) + 2*Math.PI) % (2*Math.PI);
  bendingAngle = mathRemainder(bendingAngle, 2*Math.PI);
  if (bendingAngle < 0) {
    bendingAngle = -bendingAngle;
    rotAxisAngle = ((rotAxisAngle + Math.PI) % (2*Math.PI));
  }
  if (bendingAngle < EPS || bendingAngle > Math.PI + EPS) return null;

  const numElbows    = Math.ceil(bendingAngle / maxAngle);
  const anglePerElbow = bendingAngle / numElbows;

  let allVerts = [], allFaces = [], vOffset = 0;
  let CurrentFrame = StartFrame;

  for (let i = 0; i < numElbows; i++) {
    const res = elbowGeo(r, CurrentFrame, anglePerElbow, rotAxisAngle, numSides);
    if (!res) return null;
    for (const v of res.vertices) allVerts.push(v);
    for (const [a, b, c] of res.faces)
      allFaces.push([a + vOffset, b + vOffset, c + vOffset]);
    vOffset      += res.vertices.length;
    CurrentFrame  = res.endFrame;
  }

  return { vertices: allVerts, faces: allFaces, endFrame: CurrentFrame };
}

// ── Cylinder geometry ─────────────────────────────────────────────────────────
//
// Ported from geometryHelpers.Cylinder.interpolateQtCircles

function cylinderGeo(r, start, direction, length, numSides, numCircles = 2) {
  const dhat = vunit(direction);
  const uhat = anyPerp(dhat);
  const vhat = vcross(dhat, uhat);

  const vertices = [];
  for (let i = 0; i < numCircles; i++) {
    const frac = i / (numCircles - 1);
    const p    = vadd(start, vscale(dhat, frac * length));
    for (let j = 0; j < numSides; j++) {
      const angle = 2 * Math.PI * j / numSides;
      vertices.push(vadd(p, vadd(vscale(uhat, r * Math.cos(angle)),
                                  vscale(vhat, r * Math.sin(angle)))));
    }
  }

  const faces = [];
  for (let i = 0; i < numCircles - 1; i++) {
    for (let j = 0; j < numSides; j++) {
      const nj = (j + 1) % numSides;
      const a = i*numSides + j,       b = i*numSides + nj;
      const c = (i+1)*numSides + j,   d = (i+1)*numSides + nj;
      faces.push([a, b, c]);
      faces.push([c, b, d]);
    }
  }

  return { vertices, faces };
}

// ── Main link geometry builder ────────────────────────────────────────────────

/**
 * Build CSC Dubins link geometry between two Dubins frames.
 *
 * startPose / endPose: { t:[x,y,z], x:[…], y:[…], z:[…] }
 *   (SE3 with x=col0=path-direction, y=col1, z=col2)
 *
 * prevSolution: { tUnit, tMag, c1s, c2s } for warm-start, or null.
 *
 * Returns { vertices: [[x,y,z],…], faces: [[a,b,c],…], solution } or null.
 * solution can be passed as prevSolution on the next frame.
 */
export function buildLinkGeometry(r, startPose, endPose, maxAnglePerElbow,
                                   numSides = 8, prevSolution = null) {
  const EPSILON  = 0.01;
  const DIST_EPS = r * EPSILON;

  const path = shortestCSC(r, startPose.t, startPose.x,
                                endPose.t,  endPose.x, prevSolution);

  if (vnorm(path.error) > 0.005 * r) return null;

  const solution = { tUnit: path.tUnit, tMag: path.tMag,
                     c1s: path.c1s,   c2s: path.c2s };

  let allVerts = [], allFaces = [], vOffset = 0;
  const addGeo = (verts, faces) => {
    for (const v of verts) allVerts.push(v);
    for (const [a, b, c] of faces)
      allFaces.push([a + vOffset, b + vOffset, c + vOffset]);
    vOffset += verts.length;
  };

  // ── Elbow 1  (start arc) ─────────────────────────────────────────────────
  let Elbow1EndFrame = startPose;
  if (path.theta1 > EPSILON && path.theta1 < Math.PI) {
    // rot1AxisDir = cross(sd, w1);  rot1AxisAngle = signedAngle(startPose.y, dir, sd)
    const rot1AxisDir   = vcross(startPose.x, path.w1);
    const rot1AxisAngle = signedAngle(startPose.y, rot1AxisDir, startPose.x);
    const res = compoundElbowGeo(r, startPose, path.theta1, rot1AxisAngle,
                                  maxAnglePerElbow, numSides);
    if (res) { addGeo(res.vertices, res.faces); Elbow1EndFrame = res.endFrame; }
  }

  // ── Straight cylinder ─────────────────────────────────────────────────────
  if (path.tMag > DIST_EPS) {
    const res = cylinderGeo(r, Elbow1EndFrame.t, Elbow1EndFrame.x, path.tMag, numSides);
    addGeo(res.vertices, res.faces);
  }

  // ── Elbow 2  (end arc) ───────────────────────────────────────────────────
  if (path.theta2 > EPSILON && path.theta2 < Math.PI) {
    const rot2AxisDir   = vcross(endPose.x, path.w2);
    const rot2AxisAngle = signedAngle(endPose.y, rot2AxisDir, endPose.x);

    // Elbow2StartFrame: orientation = AngleAxis(-theta2, rot2AxisDir) @ endPose.R
    //                   translation = turn2start
    const negRot = se3AngleAxis(-path.theta2, rot2AxisDir);
    const Elbow2StartFrame = {
      t: path.turn2start,
      x: mat3apply(negRot, endPose.x),
      y: mat3apply(negRot, endPose.y),
      z: mat3apply(negRot, endPose.z),
    };
    const res = compoundElbowGeo(r, Elbow2StartFrame, path.theta2, rot2AxisAngle,
                                  maxAnglePerElbow, numSides);
    if (res) addGeo(res.vertices, res.faces);
  }

  return { vertices: allVerts, faces: allFaces, solution };
}
