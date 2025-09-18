import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


poles = {
    '0 0 0 1':  ( 0,  0,  0, 1),
    '2-1-1 0':  ( 2, -1, -1, 0),
    '1-1 0 0':  ( 1, -1,  0, 0),
    '0-1 1 0':  ( 0, -1,  1, 0),
    '-1 0 1 0': (-1,  0,  1, 0),
    '-2 1 1 0': (-2,  1,  1, 0),
    '-1 1 0 0': (-1,  1,  0, 0),
    '0 1 -1 0': ( 0,  1, -1, 0),
    '1 0 -1 0': ( 1,  0, -1, 0),
    '1 1 -2 0': ( 1,  1, -2, 0),
    '-1 2 -1 0': (-1,  2, -1, 0),
    '-1 -1 2 0': (-1, -1, 2, 0),
    '1 -2 1 0': ( 1, -2,  1, 0),
}

def miller_to_miller_bravais_direction(u_p, v_p, w_p):
    u = (2 * u_p - v_p) / 3
    v = (2 * v_p - u_p) / 3
    t = -(u + v)
    w = w_p
    return [u, v, t, w]

def hkil_to_cartesian(h, k, i, l, a=1.0, c=1.633):
    """
    Convert (hkil) Miller-Bravais indices to a 3D Cartesian plane normal vector
    """
    # Basis vectors in reciprocal space (scaled)
    x = h - 0.5 * (k + i)
    y = (np.sqrt(3)/2) * (k - i)
    z = l * a / c
    return np.array([x, y, z]) / np.linalg.norm([x, y, z])


def schmid_factor(b, n, load):
    """Compute Schmid factor for a given slip system and loading direction."""
    b = b / np.linalg.norm(b)
    n = n / np.linalg.norm(n)
    load = load / np.linalg.norm(load)
    return np.abs(np.dot(load, b) * np.dot(load, n))


def miller_to_latex(vec):
    def fmt(i):
        return f"\\bar{{{abs(i)}}}" if i < 0 else str(i)
    u, v, t, w = map(int, vec)
    return f"$[{fmt(u)}{fmt(v)}{fmt(t)}{fmt(w)}]$"


def angleAxis(theta,a):
    a=a/np.linalg.norm(a) # normalize axis a
    return np.identity(3)*np.cos(theta)+np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])*np.sin(theta)+np.outer(a,a)*(1-np.cos(theta))


def normalize(v):
    v = np.array(v, dtype=float)
    norm = np.linalg.norm(v)
    if norm == 0:
        raise ValueError("Zero vector cannot be normalized.")
    return v / norm


def stereographic_projection(v):
    v = normalize(v)
    x, y, z = v
    return x / (1 + z), y / (1 + z)


def great_circle_arc(v1, v2, n=100):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    omega = np.arccos(np.clip(np.dot(v1, v2), -1, 1))
    if np.isclose(omega, 0):
        return [v1] * n
    sin_omega = np.sin(omega)
    return [(np.sin((1 - t) * omega) * v1 + np.sin(t * omega) * v2) / sin_omega for t in np.linspace(0, 1, n)]


def plot_hexagonal_poles():
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_artist(plt.Circle((0, 0), 1, edgecolor='black', facecolor='none'))
    for label, (h, k, i, l) in poles.items():
        vec = hkil_to_cartesian(h, k, i, l)
        x, y = stereographic_projection(vec)
        ax.plot(x, y, 'ko', markersize=4)
        ax.text(x+0.05, y+0.03, miller_to_latex((h, k, i, l)), fontsize=9, ha='center', va='center')
    ax.set_aspect('equal')
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.axis('off')
    plt.tight_layout()