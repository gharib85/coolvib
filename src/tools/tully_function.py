import numpy as np
from vibcooling.parser.siesta import *
#from lmfit import minimize, Parameters, Parameter, report_fit


hbar=6.58211928E-16 #eV*s
pi=3.14159265359
hbar_1=6.350780E+12  #amu*ang2/s
conversion_factor = 98226949774380.3 #for omega from sqrt(eV/amu)/Ang to Hz
me= 1822.88839 #mass of an atom relative to the mass of an electron

def tully_routine(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in eV 
#    print "     delta",delta
    gamma=0
    product=0
    hit_the_window=0
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            for i in range(len(eigenvalues[k,s,:])):
                if eigenvalues[k,s,i] > fermi_level:
                    pass
                else:
                    for f in range(len(eigenvalues[k,s,:])):
                        if eigenvalues[k,s,f] < fermi_level:
                            pass
                        elif (delta-window/2) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window/2):
                            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
                            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
#                            psi[k,s,i,:]/=np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
#                            print np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
#                            psi[k,s,f,:]/=np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
                            product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
#                                    /np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,f,:])
#                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                            gamma += np.square(np.absolute(product))*kpoints_weights[k,3]/(eigenvalues[k,s,f]-eigenvalues[k,s,i])
#                            print "Hq", H_q, "fermi_level", fermi_level, "S_q", S_q, "psi", psi[k,s,f,:]
#                            print "bracket", np.square(np.absolute(product)), "kpw", kpoints_weights[k,3]
                            hit_the_window+=1    
#                            print "gamma", gamma
    gamma *=pi*hbar_1/window#/me
    return gamma, hit_the_window

def tully_routine_eps_win(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in Hz 
#    print "     delta",delta
    gamma=0
    product=0
    window*=10
    hit_the_window=0
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
            orb_homo=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-fermi_level))
            if eigenvalues[k,s,orb_homo] > fermi_level:
                orb_homo-=1
            orb_lumo=orb_homo+1
            if eigenvalues[k,s,orb_lumo] < fermi_level:
		orb_lumo+=1
                orb_homo=orb_lumo-1
#            print fermi_level, eigenvalues[k,s,orb_homo]
            if (fermi_level-eigenvalues[k,s,orb_homo])<0:
                print "warning, homo is higher than the fermi level"
            if (fermi_level-eigenvalues[k,s,orb_lumo])>0:
                print "warning, lumo is lower than the fermi level"
            orb_min=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level-window)))
            orb_max=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level+window)))
            for i in range(orb_min,orb_lumo):
                for f in range(orb_lumo,orb_max):
                    psi[k,s,i,:]/=np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
                    psi[k,s,f,:]/=np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
                    m = np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
#                    if m != 1 + 0j:
#                    print "{:15.12f}".format(m)
#                    print np.absolute(psi[k,s,f,:]), np.absolute(psi[k,s,i,:])
#                    product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
                    product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))\
#                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))\
#                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
#                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                    gamma += np.square(np.absolute(product))*kpoints_weights[k,3]
#                    print "Hq", H_q, "fermi_level", fermi_level, "S_q", S_q, "psi", psi[k,s,f,:]
#                    print "bracket", np.square(np.absolute(product)), "kpw", kpoints_weights[k,3]
                    hit_the_window+=1     
    gamma *=pi*hbar_1/window/window#/me
    return gamma, hit_the_window


def tully_routine_eta_delta(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in Hz 
#    print "     delta",delta
    gamma=0
    product=0
    window*=1
    hit_the_window=0
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
            orb_homo=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-fermi_level))
            if eigenvalues[k,s,orb_homo] > fermi_level:
                orb_homo-=1
            orb_lumo=orb_homo+1
            if eigenvalues[k,s,orb_lumo] < fermi_level:
		orb_lumo+=1
                orb_homo=orb_lumo-1
#            print fermi_level, eigenvalues[k,s,orb_homo]
            if (fermi_level-eigenvalues[k,s,orb_homo])<0:
                print "warning, homo is higher than the fermi level"
            if (fermi_level-eigenvalues[k,s,orb_lumo])>0:
                print "warning, lumo is lower than the fermi level"
            orb_min=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level-3*window)))
            orb_max=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level+3*window)))
            for i in range(orb_min,orb_lumo):
                for f in range(orb_lumo,orb_max):
#                    if   (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (window):#delta+/-
                    if (delta-window/2) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window/2):
                        psi[k,s,i,:]/=np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
                        psi[k,s,f,:]/=np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
    #                    m = np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
    #                    if m != 1 + 0j:
    #                    print "{:15.12f}".format(m)
    #                    print np.absolute(psi[k,s,f,:]), np.absolute(psi[k,s,i,:])
    #                    product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
                        product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))\
    #                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))\
    #                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
    #                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                        gamma += np.square(np.absolute(product))*kpoints_weights[k,3]/(eigenvalues[k,s,f]-eigenvalues[k,s,i])
    #                    print "Hq", H_q, "fermi_level", fermi_level, "S_q", S_q, "psi", psi[k,s,f,:]
    #                    print "bracket", np.square(np.absolute(product)), "kpw", kpoints_weights[k,3]
                        hit_the_window+=1     
    gamma *=pi*hbar_1/window#/hit_the_window#/me
    return gamma, hit_the_window

def tully_routine_eta_fermi(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in Hz 
#    print "     delta",delta
    gamma=0
    product=0
    window*=1
    hit_the_window=0
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
            orb_homo=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-fermi_level))
            if eigenvalues[k,s,orb_homo] > fermi_level:
                orb_homo-=1
            orb_lumo=orb_homo+1
            if eigenvalues[k,s,orb_lumo] < fermi_level:
		orb_lumo+=1
                orb_homo=orb_lumo-1
#            print fermi_level, eigenvalues[k,s,orb_homo]
            if (fermi_level-eigenvalues[k,s,orb_homo])<0:
                print "warning, homo is higher than the fermi level"
            if (fermi_level-eigenvalues[k,s,orb_lumo])>0:
                print "warning, lumo is lower than the fermi level"
            orb_min=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level-3*window)))
            orb_max=min(range(len(eigenvalues[k,s,:])), key=lambda j: abs(eigenvalues[k,s,j]-(fermi_level+3*window)))
            for i in range(orb_min,orb_lumo):
                for f in range(orb_lumo,orb_max):
#                    if   (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (window):#delta+/-
                    if (delta-window/2) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window/2):
                        psi[k,s,i,:]/=np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
                        psi[k,s,f,:]/=np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
    #                    m = np.sqrt(np.dot(psi[k,s,i,:].conjugate().transpose(),psi[k,s,i,:]))
    #                    if m != 1 + 0j:
    #                    print "{:15.12f}".format(m)
    #                    print np.absolute(psi[k,s,f,:]), np.absolute(psi[k,s,i,:])
    #                    product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
                        product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))\
    #                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))\
    #                            /np.sqrt(np.dot(psi[k,s,f,:].conjugate().transpose(),psi[k,s,f,:]))
    #                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                        gamma += np.square(np.absolute(product))*kpoints_weights[k,3]/(eigenvalues[k,s,f]-eigenvalues[k,s,i])
    #                    print "Hq", H_q, "fermi_level", fermi_level, "S_q", S_q, "psi", psi[k,s,f,:]
    #                    print "bracket", np.square(np.absolute(product)), "kpw", kpoints_weights[k,3]
                        hit_the_window+=1     
    gamma *=pi*hbar_1/window#/me
    return gamma, hit_the_window




def tully_routine_k(kpoints_weights,n_spin,eigenvalues,fermi_level,omega,window,psi,H_plus,H_minus,S_plus,S_minus,dq):
    delta = hbar*omega*conversion_factor # in eV 
#    print "     delta",delta
    gamma=0
    product=0
    hit_the_window=0
    for k in range(len(kpoints_weights)):
        for s in range(n_spin):
            for i in range(len(eigenvalues[k,s,:])):
                if eigenvalues[k,s,i] > fermi_level:
                    pass
                else:
                    for f in range(len(eigenvalues[k,s,:])):
                        if eigenvalues[k,s,f] < fermi_level:
                            pass
                        elif (delta-window/2) <=  (eigenvalues[k,s,f]-eigenvalues[k,s,i]) <= (delta+window/2):
                            H_q=(H_plus[k,s,:,:]-H_minus[k,s,:,:])/2/dq
                            S_q=(S_plus[k,:,:]-S_minus[k,:,:])/2/dq
                            product= np.dot(psi[k,s,i,:].conjugate().transpose(),np.dot(H_q-fermi_level*S_q,psi[k,s,f,:]))
#                             gamma += np.square(np.absolute(product))*kpoints_weights[k,3]*np.square(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                            gamma += np.square(np.absolute(product))*kpoints_weights[k,3]/(eigenvalues[k,s,f]-eigenvalues[k,s,i])
                            hit_the_window+=1    
#                            print "gamma", gamma
    gamma *=pi*hbar_1/window
    return gamma, hit_the_window


