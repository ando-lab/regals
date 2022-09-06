#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix

class efa:
    
    def __init__(self, I, sigma = np.array([])):
        self.I = I
        self.sigma = sigma
    
    @property
    def Nr(self):
        return self.I.shape[0]
    
    @property
    def Nc(self):
        return self.I.shape[1]
    
    def evolving_factors(self, k, direction = 'forward', skip = 1):
        
        Ncols = self.Nc
        
        sv = np.nan * np.ones((k, Ncols))
        
        for j in range(0,Ncols,skip):
            if direction.lower() == 'forward':
                cols = np.arange(j+1)
            elif direction.lower() == 'reverse':
                cols = np.arange(j,Ncols)
            sv[:,j] = self.svd(k,cols)
        
        return sv
        
    def svd(self, k = None, cols = np.array([]), normalize = True, subtract = True):
        
        if k == None:
            k = self.Nc
        
        if len(cols) == 0:
            cols = np.arange(self.Nc)
        
        sv = np.nan*np.ones(k);
        
        A = self.I[:,cols]
        
        if len(self.sigma) != 0 and self.sigma.shape[0] == self.Nr:
            if len(self.sigma.shape) == 1 or self.sigma.shape[1] == 1:
                A = np.diag(1/self.sigma) @ A
            else:
                A = np.diag(1/np.mean(self.sigma[:,cols],1)) @ A
        
        if normalize:
            A = A / np.sqrt(self.Nr)
        
        s = np.linalg.svd(A, full_matrices=False, compute_uv=False)
        s = s[:min(k,len(s))]
        
        if subtract and normalize:
            s = s - np.sqrt(cols.shape[0]/self.Nr)
        elif subtract:
            s = s - np.sqrt(cols.shape[0])
        
        sv[:len(s)] = s
        return sv
    
    def quick_rotate(self, xstart, xend):
        ncomp = len(xstart)
        
        w = 1 / np.mean(self.sigma,1)
        [u, s, v] = np.linalg.svd(np.diag(w,0) @ self.I, full_matrices=False)
        
        u = u[:,:ncomp]
        s = s[:ncomp]
        v = v.T[:,:ncomp]
        
        R = np.zeros((ncomp,ncomp))
        
        for n in range(ncomp):
            m = np.full(v.shape[0],False)
            m[int(np.ceil(xstart[n])):int(np.floor(xend[n]))] = True
            
            v_in = v[m,:]
            v_out = v[~m,:]
            
            A = v_out
            B = np.mean(v_in, 0)
            AA = A.T @ A
            BB = np.outer(B,B)
            
            lmbd = 1E6 * np.trace(AA) / np.trace(BB)
            
            R[:,n] = np.linalg.solve(AA + lmbd * BB, lmbd * B * 1)
        
        y = (u @ np.diag(s) @ np.linalg.pinv(R.T)) / w[:,np.newaxis]
        c = R.T @ v.T
        
        return [y, c, R]
    
    @staticmethod
    def fit_inflection(s, direction, threshold, window):
        N = s.shape[0]
        xinfl = np.nan * np.ones(N)
        xc = np.nan *np.ones(N)
        slope = np.nan *np.ones(N)
        x = np.arange(s.shape[1])
        
        for j in range(N):
            v = s[j, :]
            is_incl = ~np.isnan(v)
            xj = x[is_incl]
            yj = v[is_incl]
            pos = np.argwhere(yj>threshold)
            if len(pos) != 0:
                if direction == 'forward':
                    x0 = xj[pos[0]]
                elif direction == 'reverse':
                    x0 = xj[pos[-1]]
                else:
                    raise ValueError('unexpected direction')
            else:
                continue
            xc[j] = x0
            is_window = np.logical_and((x0 - window / 2) <= xj, (x0 + window / 2) >= xj)
            npts = np.count_nonzero(is_window)
            if npts < 2:
                continue
            xw = xj[is_window]
            yw = yj[is_window]
            A = np.ones((npts,2))
            A[:, 0] = xw - x0
            c = np.linalg.lstsq(A, yw, rcond=None)[0]
            xinfl[j] = x0 + (1 - c[1]) / c[0]
            slope[j] = c[0]
        
        xinfl[~np.isnan(xinfl)][xinfl[~np.isnan(xinfl)]<x[0]] = x[0]
        xinfl[~np.isnan(xinfl)][xinfl[~np.isnan(xinfl)]>x[-1]] = x[-1]
        
        return [xinfl, xc, slope]
        