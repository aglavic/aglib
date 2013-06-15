#-*- coding: utf-8 -*-
'''
General purpose form factors with and without directional average.
All FF with directional average begin with capital letter.

Most function definitions are taken from:
    Rémi Lazzari, IsGISAXS: a program for graying incidence small-angle X-ray
                  scattering analysis of supported islands,
                  Journal of Applied Crystallography (2002), 35, 406-421
'''
from numpy import pi, sqrt, exp, sin, cos, sinc, ones_like
from scipy.special import j1

def Tfactor(Qx, Qy, Qz, x, y, z):
    '''
    Prefactor for any form factor to translate a particle in real space
    
    :param: x, y, z translation in real space coordinates
    '''
    return exp(1j*(Qx*x+Qy*y+Qz*z))

def sphere(Qx, Qy, Qz, R):
    '''
    Form factor of a sphere
    
    :param: R radius
    '''
    Q=sqrt(Qx**2+Qy**2+Qz**2)
    return Sphere(Q, R)

def cuboid(Qx, Qy, Qz, a, b, c):
    '''
    Form factor of a cuboid
    
    :param: a, b, c the edge lengths in x, y and z direction
    '''
    FFx=a*sinc(Qx*a/pi)
    FFy=b*sinc(Qy*b/pi)
    FFz=c*sinc(Qz*c/pi)
    return FFx*FFy*FFz

def cube(Qx, Qy, Qz, a):
    return cuboid(Qx, Qy, Qz, a, a, a)

def cylinder(Qx, Qy, Qz, R, h):
    '''
    Form factor of a cylinder oriented parallel to z-axis
    
    :param: R raidus
    :param: h height
    '''
    # xy part
    Qr=sqrt(Qx**2+Qy**2)
    QR=(Qr*R)
    FFr=ones_like(QR)
    FFrscale=2.*pi*R**2
    # make sure we don't divide by zero
    QRpos=QR!=0.
    QR=QR[QRpos]
    FFr[QRpos]=j1(QR)/QR
    # z part
    FFz=h*sinc(Qy*h/pi)
    return FFrscale*FFr*FFz

def prism(Qx, Qy, Qz, a, h):
    '''
    Form factor of a prism oriented along the z-axis and one edge parallel to x-axis
    
    :param: a edge length
    :param: h height
    '''
    # xy part
    FFxy=ones_like(Qx+Qy)
    FFxyscale=2j*(sqrt(3.))
    xypos=(Qx!=0.)&((Qx**2-Qy**2)!=0.)
    Qx=Qx[xypos]
    Qy=Qy[xypos]
    Q1=sqrt(3)*Qy-Qx
    Q2=sqrt(3)*Qy+Qx
    FFxy[xypos]=(Q1*sin(Q2*a)*exp(1j*Qx*a)-
                 Q2*sin(Q1*a)*exp(-1j*Qx*a)
                 )/(Qx*(Qx**2-Qy**2))
    # z part
    FFz=h*sinc(Qy*h/pi)
    return FFxyscale*FFxy*FFz

def truncube(Qx, Qy, Qz, a, tau):
    '''
    Form factor of a truncated cube with flat corners
    
    :param: a edge length
    :param: tau degree of truncation (0-1) 
    '''
    # Untruncated cubes form factor
    FC=cube(Qx, Qy, Qz, a)
    if tau!=0:
        a2=a/2.
        # truncated edge length is edge length/2 times tau
        b=tau*a2
        # there are a lot of cases where the corners lead to division by zero
        # so the Q-arrays are translated by a very small amount to prevent
        # either Qi=0 or Qi-Qj=0
        Qx=Qx+1e-8
        Qy=Qy+2e-8
        Qz=Qz+3e-8
        # For the truncation calculate the scattering from all 8 edges is subtracted,
        # this is done by moving and rotating a quarter of an octahedron
        # as given in By R. W. HENDRICKS, J. SCHELTEN and W. SCHMA, Philosophical Magazine (1974)
        F8=F0(Qx, Qy, Qz, b)*exp(-1j*a2*(Qx+Qy+Qz))
        F8+=F0(-Qx,-Qy,-Qz, b)*exp(-1j*a2*(-Qx-Qy-Qz))
        F8+=F0(-Qx, Qy, Qz, b)*exp(-1j*a2*(-Qx+Qy+Qz))
        F8+=F0(Qx,-Qy,-Qz, b)*exp(-1j*a2*(Qx-Qy-Qz))
        F8+=F0(Qx,-Qy, Qz, b)*exp(-1j*a2*(Qx-Qy+Qz))
        F8+=F0(-Qx, Qy,-Qz, b)*exp(-1j*a2*(-Qx+Qy-Qz))
        F8+=F0(Qx, Qy,-Qz, b)*exp(-1j*a2*(Qx+Qy-Qz))
        F8+=F0(-Qx,-Qy, Qz, b)*exp(-1j*a2*(-Qx-Qy+Qz))
        return FC-F8
    else:
        return FC

def F0(Qx, Qy, Qz, b):
    A=exp(1j*b*Qx)/(Qx*(Qx-Qy)*(Qx-Qz))
    B=exp(1j*b*Qy)/(Qy*(Qy-Qx)*(Qy-Qz))
    C=exp(1j*b*Qz)/(Qz*(Qz-Qx)*(Qz-Qy))
    D=1.0/(Qx*Qy*Qz)
    return 1j*(A+B+C-D)

############# Directional integrated variants, if available in analytic form ##############
def Sphere(Q, R):
    '''
    Form factor of a sphere
    
    :param: R radius 
    '''
    QR=Q*R
    FF=ones_like(QR)
    FFscale=4./3.*pi*R**3
    # make sure we don't divide by zero
    QRpos=QR!=0
    QR=QR[QRpos]
    FF[QRpos]=(sin(QR)-QR*cos(QR))/(QR)
    return FFscale*FF

