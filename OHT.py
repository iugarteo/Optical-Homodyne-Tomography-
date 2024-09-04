import qiskit
import numpy as np
import scipy.special as spec
import math
import matplotlib.pyplot as plt
import strawberryfields as sf
import strawberryfields.ops as ops
from strawberryfields.utils import extract_unitary, is_unitary
from itertools import product
import scipy.optimize as opt
from scipy.integrate import simps
from thewalrus.quantum import (state_vector,
                               density_matrix,
                               density_matrix_element
                               )

#__________Variables__________________________

n_cutoff=100
n_iterations = 50

name = ['2.6-squeezing-400hz', '2.76-squeezing-400hz', '2.92-squeezing-400hz', '3.14-squeezing-400hz']

#___________Functions_________________________


def _wigner_iterative(rho, l=12, g=np.sqrt(2)):
    """

    """

    x = np.linspace(-l, l, 1000)
    p = np.linspace(-l, l, 1000)
    M = np.prod(rho.shape[0])
    X, Y = np.meshgrid(x, p)
    A = 0.5 * g * (X + 1.0j * Y)

    Wlist = np.array([np.zeros(np.shape(A), dtype=complex) for k in range(M)])
    Wlist[0] = np.exp(-2.0 * np.abs(A) ** 2) / np.pi

    W = np.real(rho[0, 0]) * np.real(Wlist[0])
    for n in range(1, M):
        Wlist[n] = (2.0 * A * Wlist[n - 1]) / np.sqrt(n)
        W += 2 * np.real(rho[0, n] * Wlist[n])

    for m in range(1, M):
        temp = np.copy(Wlist[m])
        Wlist[m] = (2 * np.conj(A) * temp - np.sqrt(m) * Wlist[m - 1]) / np.sqrt(m)

        # Wlist[m] = Wigner function for |m><m|
        W += np.real(rho[m, m] * Wlist[m])

        for n in range(m + 1, M):
            temp2 = (2 * A * Wlist[n - 1] - np.sqrt(m) * temp) / np.sqrt(n)
            temp = np.copy(Wlist[n])
            Wlist[n] = temp2

            # Wlist[n] = Wigner function for |m><n|
            W += 2 * np.real(rho[m, n] * Wlist[n])

    return X, Y, 0.5 * W * g ** 2


def plot_wigner(rho, fid=None, p=None, show_log_neg=False, l=12):
    # assess spaces size

    # Q,P,W=wigner(rho,l=l)
    Q, P, W = _wigner_iterative(rho, l=l)
    boundaries = (np.min(Q), np.max(Q), np.min(P), np.max(P))
    vminmax = np.max(np.abs(W))
    plt.figure()
    # plt.grid(which='major')
    plt.imshow(np.real_if_close(W), aspect="equal", origin="lower", extent=boundaries, cmap="bwr", vmin=-vminmax,
               vmax=vminmax)
    plt.xlabel("X")
    plt.ylabel("P")
    plt.grid()
    proptext = ""
    if fid != None:
        proptext += "fidelity: %.3f \n" % fid
    if p != None:
        proptext += "probability: %.2f %% \n" % (p * 100)
    # if show_log_neg:
    # proptext+="W_log_neg: %.3f "%wigner_log_neg(rho)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    plt.gca().text(0.05, 0.95, proptext, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top',
                   bbox=props)
    plt.colorbar()
    plt.show()




"""

def fock_wave_x(x, n):
    return 1. / np.power(np.pi, 0.25) * spec.eval_hermite(n, x) / np.sqrt(2 ** n * spec.factorial(n)) * np.exp(
        -x ** 2 * 0.5)

"""

def fock_wave_x(x, N):
    fval=[]
    F_2 = math.pi**(-0.25)*np.exp(-0.5*x**2)
    F_1 = math.pi**(-0.25)*np.exp(-0.5*x**2)*x*np.sqrt(2)
    fval.append( F_2)
    fval.append( F_1)
    l=2
    while l<=N[-1]:
        F_0 = np.sqrt(2/l)*x*F_1-np.sqrt((l-1)/l)*F_2
        fval.append(F_0)
        F_2 = F_1
        F_1 = F_0
        l=l+1
    return( fval[:])


def proj_theta(Q, theta, cutoff=n_cutoff):
    fock = np.arange(cutoff)
    # vec=np.zeros(cutoff)
    vec = np.exp(1j * fock * theta) * fock_wave_x(Q, fock)
    return np.outer(vec, vec.conj())


def prob_marginal(Q, pr_theta, rho):
    return np.trace(np.dot(pr_theta, rho))


def R_iter(Q, theta, rho, probs=None):
    R = np.zeros(np.shape(rho))
    cutoff = np.shape(rho)[0]
    if probs is None:
        for q, t in zip(Q, theta):
            projector = proj_theta(q, t, cutoff=cutoff)
            pm = prob_marginal(q, projector, rho)
            if pm > 1.e-7:
                R = R + projector / prob_marginal(q, projector, rho)
    else:
        for q, t, p in zip(Q, theta, probs):
            projector = proj_theta(q, t, cutoff=cutoff)
            pm = prob_marginal(q, projector, rho)
            if pm > 1.e-7:
                R = R + p * projector / prob_marginal(q, projector, rho)
    return R


def likelihood(Q, theta, rho):
    L = 1
    for q, t in zip(Q, theta):
        L *= prob_marginal(q, proj_theta(q, t, cutoff=np.shape(rho)[0]), rho)
    return L


def G_inv(Q, theta):
    G = 0
    for q, t in zip(Q, theta):
        G = G + proj_theta(q, t)
    G_ = np.linalg.inv(G)
    return G_, G


def R_rho_R_iteration(rho, Q, theta, eps=1e20, probs=None):
    cutoff = np.shape(rho)[0]
    R = (R_iter(Q, theta, rho, probs=probs) * eps + np.eye(cutoff)) / (1 + eps)
    R_rho_R = np.dot(np.dot(R, rho), R)
    G_ = G_inv(Q, theta)[0]
    R_rho_r = np.dot(R_rho_R, G_)
    R_rho_r = np.dot(G_, R_rho_R)
    N = np.trace(R_rho_R)  # /np.power(N,1./np.shape(rho)[0])
    return R_rho_R / N



def read(n,name):
    a = open('READY_'+'%i'%n+'_%s'%name+'.txt', 'r')
    sample = []
    for line in a:
        lines = [i for i in line.split(",")]
        sample.append([float(lines[1]), float(lines[0])])
    a.close()
    return sample

#________________________________


for element in (name):
    samples = read(0,element)

    print('The data has been acquired.'+'%s'%element)

    samples = np.array(samples)
   #plt.figure()
   #plt.title("Measurement samples")
   #plt.xlabel(r"$\Theta$")
   #plt.ylabel(r"$Q_\Theta$")
   #plt.scatter(samples[:, 0], samples[:, 1]+2, s=0.5, alpha=0.5)
   #plt.show()
   #plt.savefig("Measured_state_samples %.png")


    theta = samples[:, 0]   # a phaseoffset can be added for illustration
    Q = samples[:, 1]


# do n_iterations

    rho = np.eye(n_cutoff) / float(n_cutoff)
    for i in range(n_iterations):
        print(i)
        rho = R_rho_R_iteration(rho, Q, theta)
    print('The %i iterations of the rho estimations are done.' %n_iterations)

    plt.matshow(np.abs(rho))
    plt.savefig('%s_'%element+'_reconstruction_rho.png', dpi=200)
    plt.show()
    plot_wigner(rho)
    #plt.savefig('coherent_reconstruction_wigner.png', dpi=200)

#To plot the number state population:
    number = []
    for i in range(len(rho)):
    	number.append(rho[i][i])

    cutoff = np.shape(rho)[0]
    fock = np.arange(cutoff)

    #plt.bar(fock, number, width=0.8, align='center', color='tomato')
    #plt.rcParams.update({'font.size': 12})
    #plt.grid(which='major', axis='x')
    #plt.xlabel('Photon number $n$')
    #plt.ylabel('Probability')
    #plt.savefig('SQ_ND_photon_st.png', dpi=200)
    #plt.show()

    mat = np.matrix(rho)
    string='1'+'%s'%element+'RHO.txt'
    print(string)
    a = open(string, 'w')
    for line in mat:
        np.savetxt(a, line, fmt='%.2e')
# _______________________________
print('---->DONE')