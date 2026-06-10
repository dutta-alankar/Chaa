"""Exact Riemann solver for the 1D Euler equations (Toro, ch. 4).

Used to validate the Sod shock-tube test against the exact solution.
"""
import numpy as np


def _pressure_function(p, rho, pk, ck, gamma):
    if p > pk:  # shock
        A = 2.0 / ((gamma + 1.0) * rho)
        B = (gamma - 1.0) / (gamma + 1.0) * pk
        f = (p - pk) * np.sqrt(A / (p + B))
        df = np.sqrt(A / (B + p)) * (1.0 - 0.5 * (p - pk) / (B + p))
    else:  # rarefaction
        f = 2.0 * ck / (gamma - 1.0) * ((p / pk) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        df = 1.0 / (rho * ck) * (p / pk) ** (-(gamma + 1.0) / (2.0 * gamma))
    return f, df


def star_state(rhoL, uL, pL, rhoR, uR, pR, gamma):
    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)
    p = max(0.5 * (pL + pR), 1e-8)
    for _ in range(60):
        fL, dfL = _pressure_function(p, rhoL, pL, cL, gamma)
        fR, dfR = _pressure_function(p, rhoR, pR, cR, gamma)
        dp = (fL + fR + (uR - uL)) / (dfL + dfR)
        p = max(p - dp, 1e-10)
        if abs(dp) < 1e-12 * p:
            break
    u = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    return p, u


def sample(xi, rhoL, uL, pL, rhoR, uR, pR, gamma):
    """State (rho, u, p) at similarity coordinate xi = x/t."""
    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)
    ps, us = star_state(rhoL, uL, pL, rhoR, uR, pR, gamma)
    g1 = (gamma - 1.0) / (gamma + 1.0)
    g2 = (gamma - 1.0) / (2.0 * gamma)

    if xi <= us:  # left of contact
        if ps > pL:  # left shock
            sL = uL - cL * np.sqrt((gamma + 1.0) / (2.0 * gamma) * ps / pL + g2)
            if xi <= sL:
                return rhoL, uL, pL
            rho = rhoL * ((ps / pL + g1) / (g1 * ps / pL + 1.0))
            return rho, us, ps
        # left rarefaction
        shL = uL - cL
        csL = cL * (ps / pL) ** g2
        stL = us - csL
        if xi <= shL:
            return rhoL, uL, pL
        if xi >= stL:
            return rhoL * (ps / pL) ** (1.0 / gamma), us, ps
        u = 2.0 / (gamma + 1.0) * (cL + (gamma - 1.0) / 2.0 * uL + xi)
        c = 2.0 / (gamma + 1.0) * (cL + (gamma - 1.0) / 2.0 * (uL - xi))
        rho = rhoL * (c / cL) ** (2.0 / (gamma - 1.0))
        return rho, u, rho * c * c / gamma
    # right of contact
    if ps > pR:  # right shock
        sR = uR + cR * np.sqrt((gamma + 1.0) / (2.0 * gamma) * ps / pR + g2)
        if xi >= sR:
            return rhoR, uR, pR
        rho = rhoR * ((ps / pR + g1) / (g1 * ps / pR + 1.0))
        return rho, us, ps
    # right rarefaction
    shR = uR + cR
    csR = cR * (ps / pR) ** g2
    stR = us + csR
    if xi >= shR:
        return rhoR, uR, pR
    if xi <= stR:
        return rhoR * (ps / pR) ** (1.0 / gamma), us, ps
    u = 2.0 / (gamma + 1.0) * (-cR + (gamma - 1.0) / 2.0 * uR + xi)
    c = 2.0 / (gamma + 1.0) * (cR - (gamma - 1.0) / 2.0 * (uR - xi))
    rho = rhoR * (c / cR) ** (2.0 / (gamma - 1.0))
    return rho, u, rho * c * c / gamma


def solution(x, t, x0=0.5, left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1),
             gamma=1.4):
    rho = np.empty_like(x)
    u = np.empty_like(x)
    p = np.empty_like(x)
    for n, xi in enumerate((x - x0) / t):
        rho[n], u[n], p[n] = sample(xi, *left, *right, gamma)
    return rho, u, p
