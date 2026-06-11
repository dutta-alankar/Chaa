/* Forcing.chpl — solenoidal spectral forcing with Ornstein-Uhlenbeck
 * temporal correlation (AthenaPK-style turbulence driving, simplified
 * to a modest number of discrete modes).
 *
 * The acceleration field is
 *   a(x) = sum_m cos(k_m . x + phi_m) [ A1_m e1_m + A2_m e2_m ]
 * with e1, e2 _|_ k (divergence-free), fixed random phases and OU
 * amplitudes updated once per step:
 *   A <- A exp(-dt/t_corr) + sigma sqrt(1 - exp(-2 dt/t_corr)) N(0,1).
 * Enabled with --forceAmp > 0 (Cartesian periodic boxes).
 */
module Forcing {
  use Params, Grid;
  use Math, Random;

  param maxModes = 64;
  var nModes = 0;
  var kv:  [0..#maxModes] 3*real;
  var e1v: [0..#maxModes] 3*real;
  var e2v: [0..#maxModes] 3*real;
  var phs: [0..#maxModes] real;
  var am1: [0..#maxModes] real;
  var am2: [0..#maxModes] real;

  var rng = new randomStream(real, seed = forceSeed);

  // Irwin-Hall approximation to a unit normal (good enough for driving)
  proc gauss(): real {
    var s = 0.0;
    for 1..12 do s += rng.next();
    return s - 6.0;
  }

  proc initForcing() {
    if forceAmp <= 0.0 then return;
    const L = (x1max - x1min, x2max - x2min, x3max - x3min);
    const n2hi = if act2 then forceKmax else 0;
    const n3hi = if act3 then forceKmax else 0;
    for n1 in 0..forceKmax {
      for n2 in -n2hi..n2hi {
        for n3 in -n3hi..n3hi {
          const m2 = n1*n1 + n2*n2 + n3*n3;
          if m2 < forceKmin*forceKmin || m2 > forceKmax*forceKmax then
            continue;
          // keep one of each (k, -k) pair
          if n1 == 0 && (n2 < 0 || (n2 == 0 && n3 <= 0)) then continue;
          if nModes >= maxModes then continue;

          const k = (2.0*pi*n1/L(0), 2.0*pi*n2/L(1), 2.0*pi*n3/L(2));
          const kn = sqrt(k(0)**2 + k(1)**2 + k(2)**2);
          // helper axis not parallel to k
          const h = if abs(k(2)) < 0.9*kn then (0.0, 0.0, 1.0)
                                          else (1.0, 0.0, 0.0);
          // e1 = normalize(k x h), e2 = normalize(k x e1)
          var e1 = (k(1)*h(2) - k(2)*h(1),
                    k(2)*h(0) - k(0)*h(2),
                    k(0)*h(1) - k(1)*h(0));
          const n1n = sqrt(e1(0)**2 + e1(1)**2 + e1(2)**2);
          for param c in 0..2 do e1(c) /= n1n;
          var e2 = (k(1)*e1(2) - k(2)*e1(1),
                    k(2)*e1(0) - k(0)*e1(2),
                    k(0)*e1(1) - k(1)*e1(0));
          const n2n = sqrt(e2(0)**2 + e2(1)**2 + e2(2)**2);
          for param c in 0..2 do e2(c) /= n2n;

          kv[nModes] = k;
          e1v[nModes] = e1;
          e2v[nModes] = e2;
          phs[nModes] = 2.0*pi*rng.next();
          nModes += 1;
        }
      }
    }
    if nModes == 0 then
      halt("forcing enabled but no modes selected: check forceKmin/Kmax");
    writeln("  forcing    : ", nModes, " OU modes, a_rms ~ ", forceAmp);
  }

  /* advance the OU amplitudes by dt (call once per step, serially) */
  proc updateForcing(dt: real) {
    if forceAmp <= 0.0 || nModes == 0 then return;
    const f = exp(-dt/forceTcorr);
    const sig = forceAmp/sqrt(nModes: real);
    const g = sig*sqrt(max(1.0 - f*f, 0.0));
    for m in 0..#nModes {
      am1[m] = f*am1[m] + g*gauss();
      am2[m] = f*am2[m] + g*gauss();
    }
  }

  /* acceleration at a cell centre */
  inline proc forceAccel(i: int, j: int, k: int): 3*real {
    var a: 3*real;
    const p = physPos(i, j, k);
    for m in 0..#nModes {
      const c = cos(kv[m](0)*p(0) + kv[m](1)*p(1) + kv[m](2)*p(2) + phs[m]);
      for param d in 0..2 do
        a(d) += c*(am1[m]*e1v[m](d) + am2[m]*e2v[m](d));
    }
    return a;
  }
}
