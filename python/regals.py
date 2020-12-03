#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from copy import deepcopy
import numpy as np
from scipy import sparse as sp
from scipy.linalg import eig
from scipy.sparse.linalg import spsolve

class regals:
    
    def __init__(self, I, err):
        self.I = I
        self.err = err
    
    def fit_concentrations(self, mix):
        
        mix = deepcopy(mix)
        
        H = mix.H_concentration
        [AA, Ab] = mix.concentration_problem(self.I, self.err)
        
        u = spsolve(AA + H, Ab)
        
        u = np.split(u, np.cumsum(mix.k_concentration)[:-1])
        
        mix.u_concentration = u
        return mix
    
    def fit_profiles(self, mix):
        
        mix = deepcopy(mix)
        
        H = mix.H_profile
        [AA, Ab] = mix.profile_problem(self.I, self.err)
        
        u = spsolve(AA + H, Ab)
        
        u = np.split(u, np.cumsum(mix.k_profile)[:-1])
        
        n = mix.norm_profile
        u = [uj / nj if nj != 0 else uj for uj, nj in zip(u, n)]
        
        mix.u_profile = u
        return mix
    
    def step(self, mix):
        
        new_mix = self.fit_concentrations(self.fit_profiles(mix));
        
        resid = (self.I - new_mix.I_reg) / self.err
        
        params = {}
        params['x2'] = np.mean(resid ** 2)
        params['delta_concentration'] = np.sum(np.abs(new_mix.concentrations - mix.concentrations),0)
        params['delta_profile'] = np.sum(np.abs(new_mix.profiles - mix.profiles),0)
        params['delta_u_concentration'] = np.array([np.sum(np.abs(nupk - upk)) for nupk, upk in zip(new_mix.u_concentration, mix.u_concentration)])
        params['delta_u_profile'] = np.array([np.sum(np.abs(nupr - upr)) for nupr, upr in zip(new_mix.u_profile, mix.u_profile)])
        
        return [new_mix, params, resid]
    
    def run(self, mix, stop_fun = None, update_fun = None):
        
        if stop_fun == None:
            stop_fun = lambda num_iter, params: [num_iter >= 10, 'max_iter']
        
        if update_fun == None:
            update_fun = lambda num_iter, new_mix, params, resid: True
        
        num_iter = 0
        while True:
            num_iter += 1
            [mix, params, resid] = self.step(mix)
            
            update_fun(num_iter, mix, params, resid)
            
            [if_exit, exit_cond] = stop_fun(num_iter, params)
            if if_exit:
                break
        
        return [mix, params, resid, exit_cond]



class mixture:
    
    def __init__(self, components, lambda_concentration = np.array([]), lambda_profile = np.array([]), u_concentration = [], u_profile = []):
        self.components = components
        self.lambda_concentration = lambda_concentration
        self.lambda_profile = lambda_profile
        self.u_concentration = u_concentration
        self.u_profile = u_profile
        
        if len(self.u_concentration) == 0:
            self.u_concentration = [comp.concentration.u0 for comp in components]
        
        if len(self.u_profile) == 0:
            self.u_profile = [comp.profile.u0 for comp in components]
        
        if len(self.lambda_concentration) == 0:
            self.lambda_concentration = np.zeros(self.Nc)
        
        if len(self.lambda_profile) == 0:
            self.lambda_profile = np.zeros(self.Nc)
    
    @property
    def Nc(self):
        return len(self.components)
    
    @property
    def Nx(self):
        return self.components[0].concentration.Nx
    
    @property
    def Nq(self):
        return self.components[0].profile.Nq
    
    @property
    def k_concentration(self):
        return np.array([comp.concentration.k for comp in self.components])
    
    @property
    def k_profile(self):
        return np.array([comp.profile.k for comp in self.components])
    
    @property
    def I_reg(self):
        return self.profiles @ self.concentrations.T
    
    def estimate_concentration_lambda(self,err,ng):
        AA = self.concentration_problem(np.zeros((self.Nq,self.Nx)),err,False).todense()
        
        split_pos = np.cumsum(self.k_concentration)[:-1]
        AA = [np.hsplit(AAi, split_pos) for AAi in np.vsplit(AA, split_pos)]
        
        ll = np.zeros(self.Nc)
        L = [comp.concentration.L for comp in self.components]
        for k in range(self.Nc):
            d = np.real(eig(AA[k][k], (L[k].T @ L[k]).todense())[0])
            ll[k] = mixture.ng2lambda(d, ng[k])
        
        return ll
    
    def estimate_profile_lambda(self,err,ng):
        AA = self.profile_problem(np.zeros((self.Nq,self.Nx)),err,False).todense()
        
        split_pos = np.cumsum(self.k_profile)[:-1]
        AA = [np.hsplit(AAi, split_pos) for AAi in np.vsplit(AA, split_pos)]
        
        ll = np.zeros(self.Nc)
        L = [comp.profile.L for comp in self.components]
        for k in range(self.Nc):
            d = np.real(eig(AA[k][k], (L[k].T @ L[k]).todense())[0])
            ll[k] = mixture.ng2lambda(d, ng[k])
        
        return ll
    
    def concentration_problem(self, I, err, calc_Ab = True):
        
        w = 1 / np.mean(err,1)
        
        A = [comp.concentration.A for comp in self.components]
        
        D = w[:,np.newaxis] * I
        y = self.profiles
        y = w[:,np.newaxis] * y
        
        AA = [[(y[:,k1] @ y[:,k2]) * (A[k1].T @ A[k2]) for k2 in range(self.Nc)] for k1 in range(self.Nc)]
        AA = sp.vstack((sp.hstack(tuple(AAi)) for AAi in AA))
        
        if calc_Ab == True:
            Ab = [A[k].T @ (D.T @ y[:,k]) for k in range(self.Nc)]
            Ab = np.hstack(tuple(Ab))
            return [AA, Ab]
        else:
            return AA
    
    def profile_problem(self, I, err, calc_Ab = True):
        
        w = 1 / np.mean(err,1)
        
        A = [comp.profile.A for comp in self.components]
        A = [sp.diags(w,0) @ Ai for Ai in A]
        
        D = w[:,np.newaxis] * I
        c = self.concentrations
        
        AA = [[sp.csr_matrix((c[:,k1] @ c[:,k2]) * (A[k1].T @ A[k2])) for k2 in range(self.Nc)] for k1 in range(self.Nc)]
        AA = sp.vstack((sp.hstack(tuple(AAi)) for AAi in AA))
        
        if calc_Ab == True:
            Ab = [A[k].T @ (D @ c[:,k]) for k in range(self.Nc)]
            Ab = np.hstack(tuple(Ab))
            return [AA, Ab]
        else:
            return AA
        
    def extract_concentration(self,I,err,k):
        
        notk = np.setdiff1d(np.arange(self.Nc), k)
        c = self.concentrations
        y = self.profiles
        D = I - y[:,notk] @ c[:,notk].T
        yk = y[:,k]
        w = 1/np.mean(err,1)
        m = (w**2 * yk) / (yk @ (w**2 * yk))
        pk = D.T @ m
        sigmak = np.sqrt(err.T**2 @ m**2)
        
        return [pk, sigmak]
    
    def extract_profile(self,I,err,k):
        
        notk = np.setdiff1d(np.arange(self.Nc), k)
        c = self.concentrations
        y = self.profiles
        D = I - y[:,notk] @ c[:,notk].T
        ck = c[:,k]
        m = ck/(ck @ ck)
        Ik = D @ m
        sigmak = np.sqrt(err**2 @ m**2)
        
        return [Ik, sigmak]

    @property
    def H_concentration(self):
        L = [comp.concentration.L for comp in self.components]
        B = [Lpk * (lbdpk ** 0.5) for Lpk, lbdpk in zip(L, self.lambda_concentration)]
        B = sp.block_diag(B)
        return B.T @ B
    
    @property
    def H_profile(self):
        L = [comp.profile.L for comp in self.components]
        B = [Lpr * (lbdpr ** 0.5) for Lpr, lbdpr in zip(L, self.lambda_profile)]
        B = sp.block_diag(B)
        return B.T @ B
    
    @property
    def concentrations(self):
        return np.hstack(tuple(comp.concentration.A @ upk[:,np.newaxis] for comp, upk in zip(self.components, self.u_concentration)))
    
    @property
    def profiles(self):
        return np.hstack(tuple(comp.profile.A @ upr[:,np.newaxis] for comp, upr in zip(self.components, self.u_profile)))
        
    @property
    def norm_concentration(self):
        return [comp.concentration.norm(upk) for comp, upk in zip(self.components, self.u_concentration)]
    
    @property
    def norm_profile(self):
        return [comp.profile.norm(upr) for comp, upr in zip(self.components, self.u_profile)]
    
    @staticmethod
    def ng2lambda(dd, ng):
        
        ng0 = np.count_nonzero(np.isinf(dd))
        
        dd = dd[np.logical_and(dd >= 0,~np.isinf(dd))]
        
        lambda_list = np.logspace(np.amax(np.log10(dd)) + 2, np.amin(np.log10(dd)) - 2, 51)
        
        ng_list = np.zeros(len(lambda_list))
        for j in range(len(ng_list)):
            ng_list[j] = ng0 + sum(dd / (dd + lambda_list[j]))
        
        if ng < ng_list[0]:
            optimal_lambda = np.inf
        elif ng > ng_list[-1]:
            optimal_lambda = 0
        else:
            optimal_lambda = 10 ** np.interp(ng, ng_list, np.log10(lambda_list))
        
        return optimal_lambda



class component:
    
    def __init__(self, concentration, profile):
        self.concentration = concentration
        self.profile = profile



class concentration_class:
    
    def __init__(self, concentration_type, *arg, **kwarg):
        self.concentration_type = concentration_type
        
        if self.concentration_type == 'simple':
            self._regularizer = concentration_simple(*arg, **kwarg)
        
        elif self.concentration_type == 'smooth':
            self._regularizer = concentration_smooth(*arg, **kwarg)
        
        else:
            raise ValueError('unexpected concentration type')
    
    def __getattr__(self,attr):
        return super().__getattribute__('_regularizer').__getattribute__(attr) #super() usage for deepcopy



class concentration_simple:
    
    def __init__(self, x, xmin, xmax):
        self.x = x
        self.xmin = xmin
        self.xmax = xmax
    
    @property
    def Nx(self):
        return len(self.x)
    
    @property
    def y0(self):
        return self.A @ self.u0
    
    @property
    def F(self):
        is_in_concentration = np.isin(self.x,self.w)
        ind_in_w = np.nonzero(self.x[:,np.newaxis] == self.w)[1] #assumes no repetition in x
        v = np.arange(self.Nx)
        return sp.csr_matrix((np.ones(self.Nw), (v[is_in_concentration], ind_in_w)), shape=(self.Nx, self.Nw))
    
    @property
    def k(self):
        return self.Nw
    
    @property
    def w(self):
        return self.x[np.logical_and(self.x >= self.xmin, self.x <= self.xmax)] #assumes no repetition in x
    
    @property
    def Nw(self):
        return len(self.w)
    
    @property
    def u0(self):
        return np.ones(self.k)
    
    @property
    def L(self):
        return sp.eye(self.Nw)
    
    @property
    def A(self):
        return self.F

    def norm(self,u):
        return np.mean(u**2)**0.5



class concentration_smooth:
    
    def __init__(self, x, xmin, xmax, Nw = 50, is_zero_at_xmin = True, is_zero_at_xmax = True):
        #Changed input parameter order because python requires default parameters to follow non-default ones
        self.x = x
        self.xmin = xmin
        self.xmax = xmax
        self.Nw = Nw
        self.is_zero_at_xmin = is_zero_at_xmin
        self.is_zero_at_xmax = is_zero_at_xmax
    
    @property
    def Nx(self):
        return len(self.x)
    
    @property
    def y0(self):
        return self.A @ self.u0
    
    @property
    def F(self):
        ix = np.searchsorted(self.w,self.x,side='right')
        ix = ix - 1
        is_in_concentration = np.logical_and(ix > -1, ix < self.Nw - 1)
        v = np.arange(self.Nx)
        ix = ix[is_in_concentration]
        v = v[is_in_concentration]
        u = (self.x[is_in_concentration] - self.w[ix]) * (1/self.dw)
        return sp.csr_matrix((np.concatenate((1-u,u)), (np.concatenate((v,v)), np.concatenate((ix,ix+1)))),shape = (self.Nx, self.Nw))
    
    @property
    def k(self):
        return self.Nw - self.is_zero_at_xmin - self.is_zero_at_xmax
    
    @property
    def dw(self):
        return (self.xmax - self.xmin)/(self.Nw-1)
    
    @property
    def w(self):
        return np.linspace(self.xmin, self.xmax, self.Nw)
    
    @property
    def u0(self):
        if self.is_zero_at_xmax and self.is_zero_at_xmin:
            u0_tmp =  1 - ( 2*(self.w-self.xmin)/(self.xmax-self.xmin) - 1) ** 2
            return u0_tmp[1:-1]
        elif self.is_zero_at_xmax and (not self.is_zero_at_xmin):
            u0_tmp = (self.xmax-self.w)/(self.xmax-self.xmin)
            return u0_tmp[:-1]
        elif (not self.is_zero_at_xmax) and self.is_zero_at_xmin:
            u0_tmp = (self.w-self.xmin)/(self.xmax-self.xmin)
            return u0_tmp[1:]
        else:
            return np.ones(self.Nw)
    
    @property
    def L(self):
        L_tmp = sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))
        if self.is_zero_at_xmax:
            L_tmp = L_tmp[:,:-1]
        if self.is_zero_at_xmin:
            L_tmp = L_tmp[:,1:]
        return L_tmp
    
    @property
    def A(self):
        A_tmp = self.F
        if self.is_zero_at_xmax:
            A_tmp = A_tmp[:,:-1]
        if self.is_zero_at_xmin:
            A_tmp = A_tmp[:,1:]
        return A_tmp
    
    def norm(self,u):
        return np.mean(u**2)**0.5



class profile_class:
    
    def __init__(self, profile_type, *arg, **kwarg):
        self.profile_type = profile_type
        
        if self.profile_type == 'simple':
            self._regularizer = profile_simple(*arg, **kwarg)
        
        elif self.profile_type == 'smooth':
            self._regularizer = profile_smooth(*arg, **kwarg)
        
        elif self.profile_type == 'realspace':
            self._regularizer = profile_real_space(*arg, **kwarg)
        
        else:
            raise ValueError('unexpected profile type')
    
    def __getattr__(self,attr):
        return super().__getattribute__('_regularizer').__getattribute__(attr) #super() usage for deepcopy



class profile_simple:
    
    def __init__(self, q):
        self.q = q
    
    @property
    def Nq(self):
        return len(self.q)
    
    @property
    def y0(self):
        return self.A @ self.u0
    
    @property
    def F(self):
        return sp.eye(self.Nq)
    
    @property
    def k(self):
        return self.Nq

    @property
    def u0(self):
        return np.ones(self.k)
    
    @property
    def w(self):
        return np.ones(self.k)
    
    @property
    def Nw(self):
        return self.Nq
    
    @property
    def L(self):
        return sp.eye(self.Nq)
    
    @property
    def A(self):
        return self.F
        
    def norm(self,u):
        return np.mean(u**2)**0.5



class profile_smooth:
    
    def __init__(self, q, Nw = 50):
        self.q = q
        self.Nw = Nw
    
    @property
    def Nq(self):
        return len(self.q)
    
    @property
    def y0(self):
        return self.A @ self.u0
    
    @property
    def F(self):
        ix = np.searchsorted(self.w,self.q,side='right')
        ix = ix - 1
        ix[ix == -1] = 0
        ix[ix == self.Nw - 1] = self.Nw - 2
        v = np.arange(self.Nq)
        u = (self.q - self.w[ix]) * (1/self.dw)
        return sp.csr_matrix((np.concatenate((1-u,u)), (np.concatenate((v,v)), np.concatenate((ix,ix+1)))),shape = (self.Nq, self.Nw))
    
    @property
    def k(self):
        return self.Nw
    
    @property
    def u0(self):
        return np.ones(self.k)
    
    @property
    def qmin(self):
        return np.amin(self.q)
   
    @property
    def qmax(self):
        return np.amax(self.q)
    
    @property
    def dw(self):
        return (self.qmax - self.qmin) / (self.Nw - 1)
    
    @property
    def w(self):
        return np.linspace(self.qmin, self.qmax, self.Nw)
    
    @property
    def L(self):
        return sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))
    
    @property
    def A(self):
        return self.F
    
    def norm(self,u):
        return np.mean(u**2)**0.5



class profile_real_space:
    
    def __init__(self, q, dmax, Nw = 50, is_zero_at_r0 = True, is_zero_at_dmax = True):
        self.q = q
        self.dmax = dmax
        self.Nw = Nw
        self.is_zero_at_r0 = is_zero_at_r0
        self.is_zero_at_dmax = is_zero_at_dmax
    
    @property
    def Nq(self):
        return len(self.q)
    
    @property
    def y0(self):
        return self.A @ self.u0
    
    @property
    def F(self):
        F_tmp = 4 * np.pi * self.dw * np.sinc(np.outer(self.q,self.w) / np.pi)
        F_tmp[:,[0,-1]] = 0.5*F_tmp[:,[0,-1]]
        return F_tmp
    
    @property
    def k(self):
        return self.Nw - int(self.is_zero_at_dmax) - int(self.is_zero_at_r0)
    
    @property
    def dw(self):
        return self.dmax / (self.Nw - 1)
    
    @property
    def w(self):
        return self.dmax * np.linspace(0,1,self.Nw)
    
    @property
    def u0(self):
        u0_tmp = 1 - (2 * self.w / self.dmax - 1) ** 2
        u0_tmp = u0_tmp / (4 * np.pi * self.dw * np.sum(u0_tmp))
        if self.is_zero_at_dmax:
            u0_tmp = u0_tmp[:-1]
        if self.is_zero_at_r0:
            u0_tmp = u0_tmp[1:]
        return u0_tmp
    
    @property
    def L(self):
        L_tmp = sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(0,self.Nw-2))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(1,self.Nw-1))), shape=(self.Nw-2, self.Nw))+\
            sp.csr_matrix((-0.5*np.ones(self.Nw-2), (np.arange(0,self.Nw-2), np.arange(2,self.Nw))), shape=(self.Nw-2, self.Nw))
        if self.is_zero_at_dmax:
            L_tmp = L_tmp[:,:-1]
        if self.is_zero_at_r0:
            L_tmp = L_tmp[:,1:]
        return L_tmp
    
    @property
    def A(self):
        A_tmp = self.F
        if self.is_zero_at_dmax:
            A_tmp = A_tmp[:,:-1]
        if self.is_zero_at_r0:
            A_tmp = A_tmp[:,1:]
        return A_tmp
    
    def norm(self,u):
        weight = 4 * np.pi * self.dw * np.ones(self.Nw)
        weight[[0,-1]] = 0.5* weight[[0,-1]]
        if self.is_zero_at_dmax:
            weight = weight[:-1]
        if self.is_zero_at_r0:
            weight = weight[1:]
        return np.sum(weight * u)
