! --------------------------------------------------------------------------------------------------------------
!
! Copyright (C) 2021 Pierre Mourier & Xisco Jiménez Forteza
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 3 of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

! --------------------------------------------------------------------------------------------------------------

! Computation of the (l,m,n) quasi-normal mode (QNM) complex frequency wlmn
! and angular separation constant Almn for a spin-weigth-s perturbation
! of a Kerr (or Schwarzschild) black hole using Leaver's method [1].

! This is based on the original Mathematica code from [2-4] available at [5-6];
! please do refer to these original works and note that the code available at [5-6] also provides
! a computation of the angular and radial wavefunctions of the QNMs (not included here).

! The specificities of this version are described in [7] and include:
! the use of Leaver's inversions [1] in the calculation of the
! continued fractions for a more stable recovery of overtones,
! much larger numbers of steps in the approximations to the
! continued fractions allowed by the direct implementation of the
! secant method instead of Mathematica's memory-consuming built-in root-finding algorithm,
! and a progressive increase of this number of steps over successive
! iterations along with a convergence criterion. The latter can be of use
! for modes for which Leaver's method is less efficient such as
! near the algebraically special Schwarzschild mode w = - 2j (where j = sqrt(-1)).
! The units convention also differs from [1-6] where the mass and maximal dimensionless spin are set to 1/2.
! This version requires two close but distinct initial guesses for wlmn.

! Units are set such that the mass of the black hole is 1.

! Note that this version has not been extensively tested for s /= -2 nor for l /= 2.

! --------------------------------------------------------------------------------------------------------------

! References:

! [1] E. W. Leaver, "An analytic representation for the quasi-normal modes of Kerr black holes".
! Proc. R. Soc. Lond. A 402: 285 (1985).

! [2] E. Berti, V. Cardoso and C. M. Will, "On gravitational-wave spectroscopy of massive black holes
! with the space interferometer LISA".  Phys. Rev. D 73:064030 (2006).  [arXiv: gr-qc/0512160]

! [3] E. Berti, V. Cardoso and M. Casals, "Eigenvalues and eigenfunctions of spin-weighted
! spheroidal harmonics in four and higher dimensions".  Phys. Rev. D 73:024013 (2006);
! Erratum -- ibid. 73:109902 (2006).  [arXiv: gr-qc/0511111]

! [4] E. Berti, V. Cardoso and A. O. Starinets, "TOPICAL REVIEW: Quasinormal modes of black holes
! and black branes".  Class. Quant. Grav. 26:163001 (2009).  [arXiv: 0905.2975]

! [5] V. Cardoso, https://centra.tecnico.ulisboa.pt/network/grit/files/ringdown/

! [6] E. Berti, https://pages.jh.edu/~eberti2/ringdown/

! [7] F. Jiménez Forteza and P. Mourier, "High-overtone fits to numerical relativity ringdowns:
! beyond the dismissed n=8 special tone". (2021).  [arXiv: 2107.11829]

! --------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------


program KerrQNM


implicit none

character(120) :: param_file_path

integer, parameter :: real_precision = 2   ! Set 2 to work in double precision, 4 for quadruple precision (much slower).

integer, parameter :: sp = kind(1.0), hp = selected_real_kind(real_precision*precision(1.0)), li = selected_int_kind(18)

complex(hp), parameter :: j = (0._hp, 1._hp)

integer(li) :: time1, time2
real(hp) :: timerate

integer :: l, s, m
integer :: nfreq, nang
real(hp) :: af
complex(hp) :: winit1, winit2, Almninit_alt
logical :: Almninit_default

integer(li) :: FracLevelsStart
integer :: niterMax, niterMin
real(hp) :: FracLevelsFactor, FracLevelsAngDecrFactor
real(hp) :: precisionThreshold

integer :: k1x2p1, k1pk2, ca2
real(hp) :: br

integer :: iter
integer(li) :: fraclevels, fraclevelsang
complex(hp) :: Almn
complex(hp), allocatable :: witer(:), Almniter(:)
logical :: keepconverging

real(hp) :: rootfinding_threshold


! --------------------------------------------------------------------------------------------------------------

namelist /params/ &
!
! Namelist containing the parameters of the mode, of the Kerr black hole (i.e., its spin af), and of the algorithm,
! to be set by the user in the .nml namelist file (parameter file). The path to this file is to be provided
! as a command-line argument when calling the executable.
!
! The namelist contains the following parameters:
!
            & s, &
! Spin weight of the perturbation considered (integer).
!
            & l, m, &
! l=|s|,|s|+1,|s|+2,... and m=-l,-l+1,...,l: spheroidal harmonics indices of the QNM (integers).
!
            & af, &
! Dimensionless spin of the Kerr black hole (real); 0 <= af <= 1. Setting a negative value af in [-1,0[
! is also possible and the solution obtained can be equivalently considered
! as a counter-rotating mode for a black hole of spin |af| up to wlmn --> - wlmn* and Almn --> Almn*,
! or as a corotating mode of the spheroidal harmonic mode (l,-m) for a black hole of spin |af|.
!
            & nfreq, &
! Number of inversions (integer) in the computation of Leaver's continued fraction for wlmn.
! Setting this number to the tone index n = 0,1,... of the mode looked for,
! or to a neighbouring value, should ensure a good stability of the (l,m,n) QNM solution
! for the root-finding algorithm. It does not ensure however that the solution found
! is indeed the (l,m,nfreq) mode, as the solution found is rather determined by the closest mode to the initial guesses.
!
            & nang, &
! Number of inversions (integer) in the computation of Leaver's continued fraction for Almn.
! It should be kept at nang = 0 for a more stable recovey of (l=2,m,n) modes.
! For larger l, nang > 0 may be more suitable --- e.g., nang = l-|s| seems to be appropriate in general.
!
            & winit1, winit2, &
! Two initial-guess values (complex numbers) for the complex QNM frequency, necessary for the secant method.
! They should be rather close to each other but must be distinct.
! The input format for a complex number is (a,b) where a and b are reals corresponding to the real and imaginary parts.
!
            & Almninit_default, &
! Logical (.true. or .false.). If .true., the default value Almn = l(l+1) - s(s+1) (corresponding to the Schwarzschild limit)
! will be used as the initial guess for Almn. If .false., the guess provided as Almninit_alt will be used instead.
! In most cases, there should be no need to use a value distinct from the default one.
!
            & Almninit_alt, &
! A user-defined guess for Almn (complex), to be used as its initial value if not using the default value.
! Some complex value must be provided for this parameter in any case, but it will be ignored if Almninit_default = .true..
!
            & FracLevelsStart, FracLevelsFactor, FracLevelsAngDecrFactor, niterMax, niterMin, precisionThreshold
! The continued fractions for wlmn and Almn are built at each iteration iter = 1,2,... with fraclevels(iter) terms,
! with fraclevels(iter) = FracLevelsStart*(FracLevelsFactor**(iter-1)). This value of fraclevels may optionally
! be reduced by a factor FracLevelsAngDecrFactor (e.g. setting this latter parameter to 10. or 20.)
! for the angular sector for a moderate gain of computation time, as the computation of Almn to a given
! precision typically requires much less terms in its continued fraction than the computation of wlmn.
! Once at least niterMin iterations have been performed, the loop over (wlmn, Almn) stops
! as soon as precisionThreshold is met on both the wlmn and Almn results over three consecutive iterations.
! If this convergence criterion fails to be met, the loop is stopped anyway after niterMax iterations.
! FracLevelsStart, niterMax, niterMin are integers; FracLevelsFactor, FracLevelsAngDecrFactor, precisionThreshold are reals.
!

! --------------------------------------------------------------------------------------------------------------

! Define different thresholds for considering that a zero of a function has been found, depending on whether
! double- or quadruple-precision real and complex numbers are used.
if(real_precision == 4) then
    rootfinding_threshold = 1.e-20_hp
else
    rootfinding_threshold = 1.e-13_hp
endif


call system_clock(time1)


call get_command_argument(1,param_file_path)
param_file_path = trim(param_file_path)
if(len_trim(param_file_path) == 0) then
    stop "Error: parameter file name not provided in the command-line call!"
endif
! The namelist's path is to be provided as a command-line argument, which is read by the above five lines.
! If the get_command_argument() subroutine is not recognized by the compiler (pre-Fortran 2003), comment those lines out,
! uncomment the "param_file_path = ..." line below instead and provide the namelist's path there. In that case,
! the code would need to be recompiled each time the path is changed.
!
! param_file_path = 'KerrQNMOvertones_params_example_n3.nml'


! Reading the parameters from the namelist file:
open(unit=10, file=param_file_path, status = 'old', action = 'read')
read(nml=params, unit=10)
close(10)

allocate(witer(niterMax),Almniter(niterMax))


! Defining some dependent parameters:
k1x2p1 = abs(m-s) + 1
k1pk2 = (abs(m-s) + abs(m+s)) / 2

ca2 = s*(s + 1) - k1pk2*(k1pk2 + 1)

br = sqrt(1-af**2)


! Initialise Almn to its default (Schwarzschild) value unless Almninit_default = .false., in which case use the provided value.
if (Almninit_default) then
    if(l /= s) then
        Almn = l*(l + 1) - s*(s + 1)
    else
        Almn = (0.1_hp,0.01_hp)           ! Introduce some deviation in case where Almn = 0
    endif
else
    Almn = Almninit_alt
    if(abs(Almn) < 1.D-4) then
        Almn = (0.1_hp,0.01_hp)           ! Introduce some deviation in case where Almn ~= 0
    endif
endif


! --------------------------------------------------------------------------------------------------------------

! This loop runs the successive iterations of Leaver's algorithm on the frequency and angular separation constants.
! At each iteration, the solutions for these two variables are printed successively, followed by
! the modulus of the absolute variations of the results on each variable with respect to the previous iteration.
!
iter = 1
keepconverging = .true.

do while(keepconverging)

    print '(/,A11,I6)', "Iteration:", iter
    
    fraclevels = nint(FracLevelsStart*(FracLevelsFactor**(iter-1)),li)
    fraclevelsang = nint(fraclevels/FracLevelsAngDecrFactor,li)


    witer(iter) = findrootLeaverinvw(Almn,fraclevels)
    print '(A8,(F17.12,SP,F16.12," j"))', "wlmn:  ", witer(iter)


    Almniter(iter) = findrootLeaverinva(witer(iter),Almn,fraclevelsang)
    print '(A8,(F17.12,SP,F16.12," j"))', "Almn:  ", Almniter(iter)
    
    if(iter == 1) then
        print '(A41,1P,E13.5," ;",E13.5)', "Absolute variations in wlmn and Almn:   ", &
                &  max( abs(witer(iter)-winit1), abs(witer(iter)-winit2) ), abs(Almniter(iter)-Almn)
    else 
        print '(A41,1P,E13.5," ;",E13.5)', "Absolute variations in wlmn and Almn:   ", &
                &  abs(witer(iter)-witer(iter-1)), abs(Almniter(iter)-Almniter(iter-1))
    endif
    
    Almn = Almniter(iter)


    if(iter >= niterMax) then
        keepconverging = .false.
    else if(iter >= niterMin) then
        keepconverging = ( max( abs(witer(iter) - witer(iter-1)) + abs(witer(iter-1) - witer(iter-2)), &
        & abs(Almniter(iter) - Almniter(iter-1)) + abs(Almniter(iter-1) - Almniter(iter-2)) ) >= precisionThreshold )
    endif
    
    iter = iter + 1

end do

iter = iter - 1


! --------------------------------------------------------------------------------------------------------------

! The results obtained at the successive iterations for wlmn and Almn are stored as complex values
! in witer(1), ..., witer(iter) and in Almniter(1), ..., Almniter(iter), respectively (where iter has reached its final value
!  at the end of the above loop). The last values witer(iter), Almniter(iter) are the converged results
! (unless niterMax was hit before convergence could be reached).
!
! The commands below display successively the input spin af, the solution found for wlmn, the solution found for Almn,
! the measures of convergence (total variation of the solutions obtained over the last three iterations)
! on wlmn and on Almn (separately) at the final iteration, and the total number of iterations iter at the end of the loop.
! The time taken for the completion of the calculation is then also displayed.
!
print '(//,A)', " Final results:"
if(max( abs(witer(iter) - witer(iter-1)) + abs(witer(iter-1) - witer(iter-2)), &
        & abs(Almniter(iter) - Almniter(iter-1)) + abs(Almniter(iter-1) - Almniter(iter-2)) ) >= precisionThreshold) then
    print *, "Warning: maximum number of iterations reached, but convergence is not achieved!"
endif
print '(A15,F15.12)', "af (input):   ", af
print '(A8,(F17.12,SP,F16.12," j"))', "wlmn:  ", witer(iter)
print '(A8,(F17.12,SP,F16.12," j"))', "Almn:  ", Almniter(iter)
print '(A60,1P,E13.5," ;",E13.5)', "Total variations of wlmn and of Almn over the last 3 steps:", &
                        & abs(witer(iter) - witer(iter-1)) + abs(witer(iter-1) - witer(iter-2)), &
                        & abs(Almniter(iter) - Almniter(iter-1)) + abs(Almniter(iter-1) - Almniter(iter-2))
print '(A22,I6)', "Number of iterations:", iter

call system_clock(time2, timerate)
print '(A18,F12.3,/)', "Time elapsed (s):", (time2-time1)/timerate


deallocate(witer, Almniter)


! --------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------

contains


! Computation of various coefficients appearing in Leaver's continued fractions [1].
!
function gammaw(nr,cw2,cw4)
    implicit none
    real(hp), intent(in) :: nr
    complex(hp), intent(in) :: cw2, cw4
    complex(hp) :: gammaw
    gammaw = nr**2 + (cw2 - 3)*nr + cw4 - cw2 + 2 
end function gammaw

function betaw(nr,cw1,cw3)
    implicit none
    real(hp), intent(in) :: nr
    complex(hp), intent(in) :: cw1, cw3
    complex(hp) :: betaw
    betaw = -2*(nr**2) + (cw1 + 2)*nr + cw3 
end function betaw

function alphaw(nr,cw0)
    implicit none
    real(hp), intent(in) :: nr
    complex(hp), intent(in) :: cw0
    complex(hp) :: alphaw
    alphaw = nr**2 + (cw0 + 1)*nr + cw0 
end function alphaw


function gammaa(nr,ca0)
    implicit none
    real(hp), intent(in) :: nr
    complex(hp), intent(in) :: ca0
    complex(hp) :: gammaa
    gammaa = ca0*(nr + k1pk2 + s)
end function gammaa

function betaa(nr,ca0,ca1)
    implicit none
    real(hp), intent(in) :: nr
    complex(hp), intent(in) :: ca0, ca1
    complex(hp) :: betaa
    betaa = nr*(nr-1) + 2*nr*(k1pk2 + 1 - ca0) - ca1 
end function betaa

function alphaa(nr)
    implicit none
    real(hp), intent(in) :: nr
    real(hp) :: alphaa
    alphaa = -2*(nr+1)*(nr + k1x2p1)
end function alphaa


function ccoeffswf(w,sep)
    implicit none
    complex(hp), intent(in) :: w, sep
    complex(hp) :: auxcoeff, ccoeffswf(0:4)
    auxcoeff= (2*w - af*m) / br
    ccoeffswf = (/  1 - s - 2*j*w - j*auxcoeff, -4 + 4*j*w*(2 + br) + 2*j*auxcoeff, s + 3 - 6*j*w - j*auxcoeff, &
        & (w**2)*(16 + 8*br - af**2) - 2*af*m*w - s - 1 + (4 + 2*br)*j*w - sep + (4*w + j)*auxcoeff, &
        & s + 1 - 8*(w**2) - (4*s+ 6)*j*w - (4*w + j)*auxcoeff  /)
end function ccoeffswf


function ca0f(w)
    implicit none
    complex(hp), intent(in) :: w
    complex(hp) :: ca0f
    ca0f = 2*af*w
end function ca0f
    
function ca1f(w,sep)
    implicit none
    complex(hp), intent(in) :: w, sep
    complex(hp) :: ca1f
    ca1f = 2*af*w*(k1x2p1 + s) + ca2 + (af**2)*(w**2) + sep
end function ca1f



! --------------------------------------------------------------------------------------------------------------

! Computation of Leaver's continued fraction for wlmn [1], approximated to fraclevels_cur fraction levels.
!
function Leaverfracw(w,sep,fraclevels_cur)
    implicit none
    complex(hp), intent(in) :: w, sep
    integer(li), intent(in) :: fraclevels_cur
    complex(hp) :: Leaverfracw, cwi(0:4)
    integer(li) :: n
    real(hp) :: nr
    cwi = ccoeffswf(w,sep)
    n = fraclevels_cur
    Leaverfracw = (-1._hp, 0._hp)
    do while(n > nfreq)
        nr = real(n,hp)
        Leaverfracw = gammaw(nr,cwi(2),cwi(4)) / (betaw(nr,cwi(1),cwi(3)) - alphaw(nr,cwi(0))*Leaverfracw)
        n = n - 1_li
    end do   
end function Leaverfracw

! Iteratively apply Leaver's inversion operation [1] nfreq times to the continued fraction for wlmn.
!
function Leaverinvw(w,sep,fraclevels_cur)
    implicit none
    complex(hp), intent(in) :: w, sep
    integer(li), intent(in) :: fraclevels_cur
    complex(hp) :: Leaverinvw, cwi(0:4)
    integer :: n
    real(hp) :: nr
    cwi = ccoeffswf(w,sep)
    n = 0
    Leaverinvw = cwi(3)
    do while(n < nfreq)
        nr = real(n,hp)
        Leaverinvw = betaw(nr+1,cwi(1),cwi(3)) - gammaw(nr+1,cwi(2),cwi(4))*alphaw(nr,cwi(0)) / Leaverinvw
        n = n + 1
    end do
    nr = real(nfreq, hp)
    Leaverinvw = Leaverinvw / alphaw(nr,cwi(0)) - Leaverfracw(w,sep,fraclevels_cur)
end function Leaverinvw


! Computation of Leaver's continued fraction for Almn [1], approximated to fracelevelsang_cur fraction levels.
!
function Leaverfraca(w,sep,fraclevelsang_cur)
    implicit none
    complex(hp), intent(in) :: w, sep
    integer(li), intent(in) :: fraclevelsang_cur
    complex(hp) :: Leaverfraca, ca0, ca1
    integer(li) :: n
    real(hp) :: nr
    ca0 = ca0f(w)
    ca1 = ca1f(w,sep)
    n = fraclevelsang_cur
    Leaverfraca = (-1._hp, 0._hp)
    do while(n > nang)
        nr = real(n,hp)
        Leaverfraca = gammaa(nr,ca0) / (betaa(nr,ca0,ca1) - alphaa(nr)*Leaverfraca)
        n = n - 1_li
    end do   
end function Leaverfraca

! Iteratively apply Leaver's inversion operation [1] nang times to the continued fraction for Almn.
!
function Leaverinva(w,sep,fraclevelsang_cur)
    implicit none
    complex(hp), intent(in) :: w, sep
    integer(li), intent(in) :: fraclevelsang_cur
    complex(hp) :: Leaverinva, ca0, ca1
    integer :: n
    real(hp) :: nr
    ca0 = ca0f(w)
    ca1 = ca1f(w,sep)
    n = 0
    Leaverinva = -ca1
    do while(n < nang)
        nr = real(n,hp)
        Leaverinva = betaa(nr+1,ca0,ca1) - gammaa(nr+1,ca0)*alphaa(nr) / Leaverinva
        n = n + 1
    end do
    nr = real(nang,hp)
    Leaverinva = Leaverinva / alphaa(nr) - Leaverfraca(w,sep,fraclevelsang_cur)
end function Leaverinva


! Solve for wlmn as a root of Leaver's continued fraction equation (given an estimate for Almn), using the secant method.
!
function findrootLeaverinvw(sep,fraclevels_cur)
    implicit none
    complex(hp), intent(in) :: sep
    integer(li), intent(in) :: fraclevels_cur
    complex(hp) :: findrootLeaverinvw, w1, w2, w3, fw1, fw2
    integer :: i
    i = 0
    w1 = winit1
    w2 = winit2
    fw1 = Leaverinvw(winit1,sep,fraclevels_cur)
    do while((i < 100) .and. (abs(fw1) >= rootfinding_threshold))
        fw2 = Leaverinvw(w2,sep,fraclevels_cur)
        w3 = w1 - fw1/(fw2 - fw1)*(w2 - w1)
        w1 = w2
        w2 = w3
        fw1 = fw2
        i = i+1
    end do
    findrootLeaverinvw = w1
    if(abs(fw1) >= rootfinding_threshold) then
        print *,"Warning: wlmn: root not found!"
    endif
end function findrootLeaverinvw

! Solve for Almn as a root of Leaver's continued fraction equation (given an estimate for wlmn), using the secant method.
!
function findrootLeaverinva(w,sep,fraclevelsang_cur)
    implicit none
    complex(hp), intent(in) :: w, sep
    integer(li), intent(in) :: fraclevelsang_cur
    complex(hp) :: findrootLeaverinva, Almn1, Almn2, Almn3, fAlmn1, fAlmn2
    integer :: i
    i = 0
    ! Relative variations of (+/-)1E-4 are used to generate the two initial guesses needed for Almn out of the single value
    ! obtained from the previous iteration iter-1 (or out of the single value at which Almn is initialized, if iter = 1).
    Almn1 = 0.9999_hp*sep
    Almn2 = 1.0001_hp*sep 
    fAlmn1 = Leaverinva(w,Almn1,fraclevelsang_cur)
    do while((i < 100) .and. (abs(fAlmn1) >= rootfinding_threshold))
        fAlmn2 = Leaverinva(w,Almn2,fraclevelsang_cur)
        Almn3 = Almn1 - fAlmn1/(fAlmn2 - fAlmn1)*(Almn2 - Almn1)
        Almn1 = Almn2
        Almn2 = Almn3
        fAlmn1 = fAlmn2
        i = i+1
    end do
    findrootLeaverinva = Almn1
    if(abs(fAlmn1) >= rootfinding_threshold) then
        print *,"Warning: Almn: root not found!"
    endif
end function findrootLeaverinva


! --------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------


end program KerrQNM

! --------------------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------------------
