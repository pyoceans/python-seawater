program check_values

implicit none
real*8 SP, SA, t, p, pr, longs, lats
real*8 gsw_ASal, gsw_PSal_from_ASal, gsw_alpha_t, gsw_beta_t, gsw_cp
real*8 gsw_ctmp, gsw_ptmp0_from_ctmp, gsw_dens, gsw_enthalpy, gsw_entropy
real*8  gsw_kappa, gsw_kappa_t, gsw_pden, gsw_ptmp, gsw_specvol, gsw_svel 

SP = 35.52764437773386d0; t = 25.5d0; p = 1023.d0; pr = 0.d0

SA = 35.7d0

longs = 201.d0; lats = -21.d0

print *
print *, 'INPUT'
print *, '*****'; print * 
print *, 'SP    = ', SP; print * 
print *, 'SA    = ', SA; print * 
print *, 't     = ', t;  print *
print *, 'p     = ', p;  print *
print *, 'pr    = ', pr; print *
print *, 'longs = ', longs; print *
print *, 'lats  = ', lats

!******************************************************************************
    
print *; print *; print *; print *, 'OUTPUT'
print *, '******'

print *, '                            computed value          expected computed value'
    
!******************************************************************************

print*; print *, 'gsw_ASal               : ', gsw_ASal(SP,p,longs,lats), 35.7d0

print*; print *, 'gsw_PSal_from_ASal     : ', gsw_PSal_from_ASal(SA,p,longs,lats), 35.5276443777339d0

print*; print *, 'gsw_alpha_t            : ', gsw_alpha_t(SA,t,p), 0.0003098378393192645d0

print*; print *, 'gsw_beta_t             : ', gsw_beta_t(SA,t,p), 0.0007257297978386655d0

print*; print *, 'gsw_cp                 : ', gsw_cp(SA,t,p), 3974.42541259729d0
           
print*; print *, 'gsw_ctmp               : ', gsw_ctmp(SA,t), 25.4805463842239d0 
           
print*; print *, 'gsw_ptmp0_from_ctmp    : ', gsw_ptmp0_from_ctmp(SA,25.4805463842239d0), 25.5d0 

print*; print *, 'gsw_density            : ', gsw_dens(SA,t,p),1027.95249315662d0
           
print*; print *, 'gsw_enthalpy           : ', gsw_enthalpy(SA,t,p), 110776.712408975d0
                      
print*; print *, 'gsw_entropy            : ', gsw_entropy(SA,t,p), 352.81879771528d0
                      
print*; print *, 'gsw_kappa              : ', gsw_kappa(SA,t,p), 4.033862685464779d-6
                                 
print*; print *, 'gsw_kappa_t            : ', gsw_kappa_t(SA,t,p), 4.104037946151349d-6
                                 
print*; print *, 'gsw_pden               : ', gsw_pden(SA,t,p,pr), 1023.66254941185d0       
                                 
print*; print *, 'gsw_ptmp               : ', gsw_ptmp(SA,t,p,pr), 25.2720983155409d0       
                                 
print*; print *, 'ptmp_inverse           : ', gsw_ptmp(SA,25.2720983155409d0,pr,p), 25.5d0 

print*; print *, 'gsw_specvol            : ', gsw_specvol(SA,t,p), 0.0009728076021579713d0 
           
print*; print *, 'gsw_svel               : ', gsw_svel(SA,t,p), 1552.93372863425d0 


print *; print *; pause

stop
end

           