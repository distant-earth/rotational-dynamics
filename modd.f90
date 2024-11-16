module modd
implicit none

! DOP853 intrinsic parameters:
integer :: N
integer :: LWORK
integer :: LIWORK
real(8), allocatable :: WORK(:)
integer, allocatable :: IWORK(:)
real(8), allocatable :: y(:)
real(8) :: RTOL
real(8) :: ATOL
real(8) :: ITOL
real(8) :: RPAR
integer :: IPAR
integer :: IOUT
integer :: IDID

! Constants:
real(8), parameter :: GM = 19.909885359645706d0 ! [earthRad**3 / h**2]
real(8), parameter :: pi = 4 * atan(1.0d0)

! Problem parameters set by the user:
real(8) :: P_0 ! rotational period [h]
real(8) :: q, ecc ! orbital parameters: pericentric distance [earthRad] and eccentricity
real(8) :: phi_0, theta_0, psi_0 ! initial Euler angles [rad]
real(8) :: AC, BC ! moments of inertia ratios
real(8) :: sphere ! numerical integration is performed within a sphere of this radius [earthRad]

! Problem parameters to be calculated:
real(8) :: a ! semimajor axis of hyperbolic orbit, [au]
real(8) :: n_mean ! mean motion, [1 / h]
real(8) :: omega_0(3) ! rotational velocity vector, [1 / h]
real(8) :: H0 ! initial mean anomaly of an asteroid on a hyperbolic orbit, [rad]
real(8) :: time ! integration time within a sphere of 'sphere' radii, [h]


contains

subroutine set_dop853_params()
	N = 7 ! number of equations
	LWORK = 11 * N + 20
	LIWORK = 20
	allocate(WORK(LWORK), IWORK(LIWORK), y(N))
	! Local error is expected to be less than RTOL*ABS(Y(I))+ATOL if ITOL = 0
	! and less than RTOL(I)*ABS(Y(I))+ATOL(I) if ITOL = 1:
	RTOL = 1e-14
	ATOL = 1e-14
	! Set ITOL = 0 if RTOL and ATOL are scalars,
	! set ITOL = 1 if RTOL and ATOL are vectors:
	ITOL = 0
	! Set IOUT = 1 if 'solout' subroutine provides output,
	! set IOUT = 0 if no output is needed:
	IOUT = 1
	! Set working parameters to default (0):
	WORK = 0.d0
	IWORK = 0
	! Specify working parameters (if needed):
	IWORK(1) = 2147483647 ! max number of steps (recommended)
	!WORK(3) = 0.333d0 ! upper limit for new_step / old_step ratio
	!WORK(4) = 6.d0 ! lower limit for new_step / old_step ratio
	!WORK(6) = 1e-4 ! max step size (if set to 0 then max step size is t_final - t_0)
	!WORK(7) = 0.d0 ! initial step size
end subroutine set_dop853_params


subroutine init_params()
	open(unit = 0, file = 'INPUT')
	read(0,*) P_0
	read(0,*) q
	read(0,*) ecc
	read(0,*) phi_0 
	phi_0 = phi_0 * pi / 180.d0
	read(0,*) theta_0
	theta_0 = theta_0 * pi / 180.d0
	read(0,*) psi_0
	psi_0 = psi_0 * pi / 180.d0
	read(0,*) AC
	read(0,*) BC
	read(0,*) sphere
end subroutine init_params


subroutine set_params()
real(8) :: cos_f, cos_f2, ksi
	a = q / (ecc - 1)
	n_mean = sqrt(GM) * a**(-1.5)
	omega_0(1) = 0.d0
	omega_0(2) = 0.d0
	omega_0(3) = 2 * pi / P_0
	cos_f = (a * (ecc**2 - 1) / sphere - 1) / ecc
	cos_f2 = (1 + cos_f) / 2
	ksi = - sqrt((ecc - 1) / (ecc + 1)) * sqrt(1 / cos_f2 - 1)
	H0 = log((1 + ksi) / (1 - ksi))
	time = 2 * (H0 - ecc * sinh(H0)) / n_mean
end subroutine set_params


function Euler2RH(phi, theta, psi)
real(8) :: phi, theta, psi
real(8) :: Euler2RH(4)
	Euler2RH(1) = cos(theta / 2) * cos(phi / 2) * cos(psi / 2) - sin(theta / 2) * sin(phi / 2) * sin(psi / 2)
	Euler2RH(2) = cos(theta / 2) * sin(phi / 2) * cos(psi / 2) - sin(theta / 2) * cos(phi / 2) * sin(psi / 2)
	Euler2RH(3) = cos(theta / 2) * cos(phi / 2) * sin(psi / 2) + sin(theta / 2) * sin(phi / 2) * cos(psi / 2)
	Euler2RH(4) = cos(theta / 2) * sin(phi / 2) * sin(psi / 2) + sin(theta / 2) * cos(phi / 2) * cos(psi / 2)
	return
end function Euler2RH


function RH2Euler(l0, l1, l2, l3)
real(8) :: l0, l1, l2, l3, RH2Euler(3)
	RH2Euler(1) = asin(2 * (l0 * l1 + l2 * l3)) ! phi
	RH2Euler(2) = atan2(2 * (l0 * l3 - l1 * l2), 1 - 2 * (l1**2 + l3**2)) ! theta
	RH2Euler(3) = atan2(2 * (l0 * l2 - l1 * l3), 1 - 2 * (l1**2 + l2**2)) ! psi
end function RH2Euler


function coeff(AC, BC)
real(8) :: AC, BC, AB
real(8) :: coeff(3)
	AB = AC / BC
	coeff(1) = 1 / AB - 1 / AC
	coeff(2) = 1 / BC - AB
	coeff(3) = AC - BC
	return
end function coeff


function abc2xyz(vec_abc, l0, l1, l2, l3)
real(8) :: vec_abc(3), abc2xyz(3)
real(8) :: R_a(3,3), R_b(3,3), R_c(3,3)
real(8) :: l0, l1, l2, l3
real(8) :: euler(3), phi, theta, psi
	euler = RH2Euler(l0, l1, l2, l3)
	phi = euler(1)
	theta = euler(2)
	psi = euler(3)
	R_a = transpose(reshape((/ 1.d0, 0.d0, 0.d0, 0.d0, cos(phi), -sin(phi), 0.d0, sin(phi), cos(phi) /), shape(R_a)))
	R_b = transpose(reshape((/ cos(psi), 0.d0, sin(psi), 0.d0, 1.d0, 0.d0, -sin(psi), 0.d0, cos(psi) /), shape(R_b)))
	R_c = transpose(reshape((/ cos(theta), -sin(theta), 0.d0, sin(theta), cos(theta), 0.d0, 0.d0, 0.d0, 1.d0 /), shape(R_c)))
	abc2xyz = matmul(R_c, matmul(R_a, matmul(R_b, vec_abc)))
	return
end function abc2xyz


! Right parts of equations:
subroutine func_f(N, t, y, res, RPAR, IPAR)
integer :: N, IPAR
real(8) :: RPAR
real(8) :: res(N), y(N)
real(8) :: t
real(8) :: r, f
real(8) :: alpha, beta, ggamma
real(8) :: omega_a, omega_b, omega_c
real(8) :: lam0, lam1, lam2, lam3
real(8) :: cos_ff, sin_ff
real(8) :: H, H_old, M, delta
real(8) :: coefs(3)
	omega_a = y(1)
	omega_b = y(2)
	omega_c = y(3)
	lam0 = y(4)
	lam1 = y(5)
	lam2 = y(6)
	lam3 = y(7)
	! Solving Kepler equation:
	M = n_mean * t
	H_old = 0.d0
	delta = 1
	do while (delta >= 1e-14)
		H = asinh((H_old + M) / ecc)
		delta = abs(H - H_old)
		H_old = H
	enddo
	f = 2 * atan(sqrt((ecc + 1) / (ecc - 1)) * tanh(H / 2))
	alpha = (lam0**2 + lam1**2 - lam2**2 - lam3**2) * cos(f)+ 2 * (lam0 * lam3 + lam1 * lam2) * sin(f)
	beta = 2 * (lam1 * lam2 - lam0 * lam3) * cos(f) + (lam0**2 - lam1**2 + lam2**2 - lam3**2) * sin(f)
	ggamma = 2 * ((lam0 * lam2 + lam1 * lam3) * cos(f) + (lam2 * lam3 - lam0 * lam1) * sin(f))
	r = a * (ecc**2 - 1) / (1 + ecc * cos(f))
	coefs = coeff(AC, BC)
	res(1) = coefs(1) * (omega_b * omega_c - 3 * GM * beta * ggamma / r**3)
	res(2) = coefs(2) * (omega_c * omega_a - 3 * GM * alpha * ggamma / r**3)
	res(3) = coefs(3) * (omega_b * omega_a - 3 * GM * beta * alpha / r**3)
	res(4) = - 0.5 * (lam1 * omega_a + lam2 * omega_b + lam3 * omega_c)
	res(5) = 0.5 * (lam0 * omega_a - lam3 * omega_b + lam2 * omega_c)
	res(6) = 0.5 * (lam3 * omega_a + lam0 * omega_b - lam1 * omega_c)
	res(7) = - 0.5 * (lam2 * omega_a - lam1 * omega_b - lam0 * omega_c)
end subroutine func_f


! Subroutine for output:
subroutine solout(NR, t_old, t, y, N, con, icomp, ND, RPAR, IPAR, IRTRN)
implicit none
integer :: N, NR, ND, IPAR, IRTRN
real(8) :: y(N), t_old, t, con, icomp, RPAR
real(8) :: w1, w2, w3, l0, l1, l2, l3
real(8) :: w_abc(3), w_xyz(3), ort_z(3)
real(8) :: new_angle
	w1 = y(1)
	w2 = y(2)
	w3 = y(3)
	l0 = y(4)
	l1 = y(5)
	l2 = y(6)
	l3 = y(7)
	w_abc = [w1, w2, w3]
	w_xyz = abc2xyz(w_abc, l0, l1, l2, l3)
	ort_z = [0.d0, 0.d0, 1.d0]
	new_angle = acos(dot_product(ort_z, w_xyz) / norm2(w_xyz))
	write(1, *) t, 2 * pi / norm2(w_abc), new_angle, 'norm:', norm2([l0, l1, l2, l3])
end subroutine solout

end module modd
