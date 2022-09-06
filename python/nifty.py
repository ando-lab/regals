#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import sparse as sp

class nifty:
    
    def __init__(self):
        self.q = None
        self.intensity = None
        self.err = None
        self.dmax = None
        self.Nr = 50
        self.is_zero_at_dmax = True
        self.is_zero_at_zero = True
        
    @property
    def F(self):
        F_tmp = 4 * np.pi * self.dr * np.sinc(np.outer(self.q,self.r) / np.pi)
        F_tmp[:,[0,-1]] = 0.5*F_tmp[:,[0,-1]]
        return F_tmp
    
    @property
    def k(self):
        return self.Nr - int(self.is_zero_at_dmax) - int(self.is_zero_at_zero)
    
    @property
    def dr(self):
        return self.dmax / (self.Nr - 1)
    
    @property
    def r(self):
        return self.dmax * np.linspace(0,1,self.Nr)
    
    @property
    def Nq(self):
        return len(self.q)
    
    @property
    def L(self):
        L_tmp = sp.csr_matrix((-0.5*np.ones(self.Nr-2), (np.arange(0,self.Nr-2), np.arange(0,self.Nr-2))), shape=(self.Nr-2, self.Nr))+\
            sp.csr_matrix((np.ones(self.Nr-2), (np.arange(0,self.Nr-2), np.arange(1,self.Nr-1))), shape=(self.Nr-2, self.Nr))+\
            sp.csr_matrix((-0.5*np.ones(self.Nr-2), (np.arange(0,self.Nr-2), np.arange(2,self.Nr))), shape=(self.Nr-2, self.Nr))
        if self.is_zero_at_dmax:
            L_tmp = L_tmp[:,:-1]
        if self.is_zero_at_zero:
            L_tmp = L_tmp[:,1:]
        return L_tmp
    
    @property
    def A(self): #Did not consider the case when err is a constant
        A_tmp = self.F
        if len(self.err) != 0:
            A_tmp = (1 / self.err[:,np.newaxis] * np.ones(self.Nr)) * A_tmp
        if self.is_zero_at_dmax:
            A_tmp = A_tmp[:,:-1]
        if self.is_zero_at_zero:
            A_tmp = A_tmp[:,1:]
        return A_tmp
    
    @property
    def b(self):
        return self.intensity/self.err
    
    def ift(self, alpha=None):
        Lmat = self.L
        Amat = self.A
        AA = Amat.T @ Amat
        LL = sp.csr_matrix.todense(Lmat.T @ Lmat)
        if alpha == None:
            alpha = np.trace(AA)/np.trace(LL)
        w = np.linalg.lstsq(AA + alpha * LL, Amat.T @ self.b, rcond=None)[0]
        if self.is_zero_at_dmax:
            w = np.append(w, 0)
        if self.is_zero_at_zero:
            w = np.insert(w, 0, 0)
        return [w, alpha]
    
    @staticmethod
    def sprite(q,intensity,err,dmax,Nr = 50,is_zero_at_dmax = True,is_zero_at_zero = True):
        N = nifty()
        N.q = q
        N.intensity = intensity
        N.err = err
        N.dmax = dmax
        N.Nr = Nr
        N.is_zero_at_dmax = is_zero_at_dmax
        N.is_zero_at_zero = is_zero_at_zero
        
        [w, alpha] = N.ift()
        r = N.r
        ireg = N.F @ w
        
        return [r, w, ireg, alpha]
        
