/* =========================================================
   shared-math-tools.js
   Shared math helpers used by mat225.html / mat227.html / mat276.html
========================================================= */

/* =========================
   DOM helpers
========================= */
export const $ = (id) => document.getElementById(id);
export function isFiniteAll(arr){ return arr.every(Number.isFinite); }

/* =========================
   Complex arithmetic
========================= */
export function complex(re, im=0){ return {re, im}; }
export function cAdd(x,y){ return complex(x.re+y.re, x.im+y.im); }
export function cSub(x,y){ return complex(x.re-y.re, x.im-y.im); }
export function cMul(x,y){ return complex(x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re); }
export function cDiv(x,y){
  const den = y.re*y.re + y.im*y.im;
  return complex((x.re*y.re + x.im*y.im)/den, (x.im*y.re - x.re*y.im)/den);
}
export function cAbs(x){ return Math.hypot(x.re, x.im); }
export function cIsZero(x, eps=1e-12){ return cAbs(x) < eps; }

export function fmtComplex(z, digits=6){
  const re = Math.abs(z.re) < 1e-12 ? 0 : z.re;
  const im = Math.abs(z.im) < 1e-12 ? 0 : z.im;
  if (im === 0) return `${re.toFixed(digits).replace(/\.?0+$/,"")}`;
  const sign = im >= 0 ? "+" : "-";
  return `${re.toFixed(digits).replace(/\.?0+$/,"")} ${sign} ${Math.abs(im).toFixed(digits).replace(/\.?0+$/,"")}i`;
}

export function fmtVecComplex(v){
  const allReal = v.every(z => Math.abs(z.im) < 1e-12);
  if (allReal) {
    const xs = v.map(z => (Math.abs(z.re) < 1e-12 ? 0 : z.re).toFixed(6));
    return `(${xs.join(", ")})`;
  }
  let maxIdx = 0;
  for (let i=1;i<v.length;i++) if (cAbs(v[i]) > cAbs(v[maxIdx])) maxIdx = i;
  const scale = cIsZero(v[maxIdx]) ? complex(1,0) : v[maxIdx];
  const vv = v.map(z => cDiv(z, scale));
  return `(${vv.map(z => fmtComplex(z)).join(", ")})`;
}

/* =========================
   2×2 eigen helper
========================= */
export function eigen2x2(a,b,c,d){
  const tr = a + d;
  const det = a*d - b*c;
  const disc = tr*tr - 4*det;

  const sqrtDisc = disc >= 0 ? complex(Math.sqrt(disc), 0) : complex(0, Math.sqrt(-disc));
  const halfTr = complex(tr/2, 0);
  const halfS  = cMul(sqrtDisc, complex(0.5,0));

  const lam1 = cAdd(halfTr, halfS);
  const lam2 = cSub(halfTr, halfS);

  function eigenvector2x2(lam){
    const aa = cSub(complex(a,0), lam);
    const dd = cSub(complex(d,0), lam);
    const bb = complex(b,0);
    const cc = complex(c,0);

    const r1 = cAbs(aa) + cAbs(bb);
    const r2 = cAbs(cc) + cAbs(dd);

    let x,y;
    if (r1 >= r2) { x = cMul(complex(-1,0), bb); y = aa; }
    else { x = cMul(complex(-1,0), dd); y = cc; }

    if (cIsZero(x) && cIsZero(y)) return [complex(1,0), complex(0,0)];
    return [x,y];
  }

  return { lam1, lam2, v1: eigenvector2x2(lam1), v2: eigenvector2x2(lam2) };
}

/* =========================
   3×3 eigen (Cardano cubic + row-cross eigenvector)
========================= */
export function det3(A){
  const [[a,b,c],[d,e,f],[g,h,i]] = A;
  return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
}
export function mul3(A,B){
  const C = Array.from({length:3}, ()=>[0,0,0]);
  for (let r=0;r<3;r++) for (let c=0;c<3;c++){
    let s = 0;
    for (let k=0;k<3;k++) s += A[r][k]*B[k][c];
    C[r][c]=s;
  }
  return C;
}
export function tr3(A){ return A[0][0] + A[1][1] + A[2][2]; }
export function trA2(A){ return tr3(mul3(A,A)); }
export function realCbrt(x){ return x >= 0 ? Math.pow(x, 1/3) : -Math.pow(-x, 1/3); }

export function cubicRootsRealCoeffs(t, s, d){
  // λ^3 - t λ^2 + s λ - d = 0  (real coefficients)
  const a = -t, b = s, c = -d;
  const p = b - a*a/3;
  const q = 2*a*a*a/27 - a*b/3 + c;
  const halfQ = q/2;
  const thirdP = p/3;
  const Delta = halfQ*halfQ + thirdP*thirdP*thirdP;

  if (Delta >= 0) {
    const sqrtD = Math.sqrt(Delta);
    const u = realCbrt(-halfQ + sqrtD);
    const v = realCbrt(-halfQ - sqrtD);
    const u1 = u+v;

    const re2 = -(u+v)/2;
    const im2 = (Math.sqrt(3)/2)*(u - v);

    return [
      complex(u1 - a/3, 0),
      complex(re2 - a/3, im2),
      complex(re2 - a/3, -im2)
    ];
  } else {
    const r = Math.sqrt(-thirdP);
    const phi = Math.acos((-halfQ)/(r*r*r));
    const x1 = 2*r*Math.cos(phi/3) - a/3;
    const x2 = 2*r*Math.cos((phi + 2*Math.PI)/3) - a/3;
    const x3 = 2*r*Math.cos((phi + 4*Math.PI)/3) - a/3;
    return [complex(x1,0), complex(x2,0), complex(x3,0)];
  }
}

export function rowCrossEigenvector3(A, lam){
  const M = [
    [cSub(complex(A[0][0],0), lam), complex(A[0][1],0),            complex(A[0][2],0)],
    [complex(A[1][0],0),            cSub(complex(A[1][1],0), lam), complex(A[1][2],0)],
    [complex(A[2][0],0),            complex(A[2][1],0),            cSub(complex(A[2][2],0), lam)]
  ];
  function rowNorm(r){ return cAbs(r[0]) + cAbs(r[1]) + cAbs(r[2]); }
  const idx = [0,1,2].sort((i,j)=>rowNorm(M[j]) - rowNorm(M[i]));
  const pairs = [[idx[0],idx[1]],[idx[0],idx[2]],[idx[1],idx[2]]];

  function cross(u,v){
    return [
      cSub(cMul(u[1], v[2]), cMul(u[2], v[1])),
      cSub(cMul(u[2], v[0]), cMul(u[0], v[2])),
      cSub(cMul(u[0], v[1]), cMul(u[1], v[0]))
    ];
  }
  for (const [i,j] of pairs){
    const w = cross(M[i], M[j]);
    if (!cIsZero(w[0],1e-10) || !cIsZero(w[1],1e-10) || !cIsZero(w[2],1e-10)) return w;
  }
  return [complex(1,0), complex(0,0), complex(0,0)];
}

/* =========================
   2×2 linear system solver
========================= */
export function linear2_solve(prefix){
  const a1 = parseFloat($(prefix+"_a1").value);
  const b1 = parseFloat($(prefix+"_b1").value);
  const c1 = parseFloat($(prefix+"_c1").value);
  const a2 = parseFloat($(prefix+"_a2").value);
  const b2 = parseFloat($(prefix+"_b2").value);
  const c2 = parseFloat($(prefix+"_c2").value);

  if (!isFiniteAll([a1,b1,c1,a2,b2,c2])) return "Invalid input.";

  const det = a1*b2 - a2*b1;
  const eps = 1e-12;

  if (Math.abs(det) > eps) {
    const x = (c1*b2 - c2*b1) / det;
    const y = (a1*c2 - a2*c1) / det;
    return `x = ${x.toFixed(10).replace(/\.?0+$/,"")}\ny = ${y.toFixed(10).replace(/\.?0+$/,"")}`;
  }

  function isZero(z){ return Math.abs(z) < eps; }
  const r1 = [a1,b1,c1], r2 = [a2,b2,c2];

  let k = null;
  for (let i=0;i<3;i++){ if (!isZero(r1[i])) { k = r2[i]/r1[i]; break; } }

  if (k === null) {
    return isZero(c1) ? "Infinitely many solutions." : "No solution.";
  }

  const proportional =
    Math.abs(r2[0] - k*r1[0]) < 1e-9 &&
    Math.abs(r2[1] - k*r1[1]) < 1e-9 &&
    Math.abs(r2[2] - k*r1[2]) < 1e-9;

  return proportional ? "Infinitely many solutions." : "No solution.";
}

/* =========================
   MAT276 equilibrium (numerical)
========================= */
export function makeSafeFunction(expr){
  const body = `
    "use strict";
    const {
      abs, acos, asin, atan, atan2, ceil, cos, exp, floor, log,
      max, min, pow, round, sin, sqrt, tan, PI, E
    } = Math;
    return (${expr});
  `;
  return new Function("x","y", body);
}

export function newton2D(f, g, x0, y0, maxIters, tol){
  let x = x0, y = y0;
  const h = 1e-6;

  for (let k=0; k<maxIters; k++){
    const F = f(x,y), G = g(x,y);
    const r = Math.hypot(F,G);
    if (!Number.isFinite(r)) return null;
    if (r < tol) return {x, y, r, iters:k};

    const Fx = (f(x+h, y) - f(x-h, y)) / (2*h);
    const Fy = (f(x, y+h) - f(x, y-h)) / (2*h);
    const Gx = (g(x+h, y) - g(x-h, y)) / (2*h);
    const Gy = (g(x, y+h) - g(x, y-h)) / (2*h);

    const detJ = Fx*Gy - Fy*Gx;
    if (!Number.isFinite(detJ) || Math.abs(detJ) < 1e-14) return null;

    const dx = (-F*Gy + Fy*G) / detJ;
    const dy = ( -Fx*G + F*Gx) / detJ;

    const step = Math.hypot(dx,dy);
    const alpha = step > 1 ? 1/step : 1;

    x += alpha*dx;
    y += alpha*dy;
    if (!Number.isFinite(x) || !Number.isFinite(y)) return null;
  }

  const rFinal = Math.hypot(f(x,y), g(x,y));
  if (Number.isFinite(rFinal) && rFinal < tol) return {x, y, r:rFinal, iters:maxIters};
  return null;
}

export function dedupPush(roots, cand, eps){
  for (const r of roots){
    if (Math.hypot(r.x - cand.x, r.y - cand.y) < eps) return false;
  }
  roots.push(cand);
  return true;
}
