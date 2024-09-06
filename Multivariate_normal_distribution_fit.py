import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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


def plot_wigner(rho, fid=None, p=None, show_log_neg=False, l=6):
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


#reconstructed state to be loaded:
name='5sqisolator800hzRHO.txt'

a= open('%s'%name, 'r')
A=np.loadtxt(a, complex)


#plot_wigner(A)

X,Y,W=_wigner_iterative(A)


X=X.ravel()
Y=Y.ravel()
W=W.ravel()
yy=[X,Y]
data_ref=np.array(yy)

#If a p0 is not given the default values are 1 and it raises an error since 1/infty will appear in the function
p0=[0.5, 0.5, 0.1, 0, 0]

def multivariate(yy,sigmax, sigmay,p, mux, muy):
    x,y=yy
    frac=1./(2*np.pi*sigmax*sigmay*np.sqrt(1-p**2))
    other=-1./(2*(1-p**2))*(((x-mux)/sigmax)**2+((y-muy)/sigmay)**2-2*p*((x-mux)/sigmax)*((y-muy)/sigmay))
    return frac*np.exp(other)


fit = curve_fit(multivariate, data_ref,W, p0)
sigmax=fit[0][0]
sigmay=fit[0][1]
p=fit[0][2]
mux=fit[0][3]
muy=fit[0][4]
mean=[mux, muy]

print('The mean vector: ', mean)
covariance=[[sigmax**2, p*sigmax*sigmay], [p*sigmax*sigmay, sigmay**2]]
print('The covariance matrix: ', covariance)

#The matrix is diagonalized to obtain the squeezing and antisqueezing values
eigenvalues, eigenvectors=np.linalg.eig(covariance)
D=np.diag(eigenvalues)
print('$\sigma^2_x(W)=$', eigenvalues[0])
print('$\sigma^2_y(W)=$', eigenvalues[1])