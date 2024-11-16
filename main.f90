program main
use modd
implicit none
real(8) :: t_0, t_final
real(8) :: euler(3)
	call set_dop853_params()
	call init_params()
	call set_params()
	t_0 = - time / 2
	t_final = time / 2
	y(1:3) = omega_0
	y(4:7) = Euler2RH(phi_0, theta_0, psi_0)
	open(unit = 1, file = 'RESULT')
	write(1, *) ecc, q, phi_0, P_0, AC, BC, '# ecc, q, phi_0, P_0, A/C, B/C'
	call DOP853(N, func_f, t_0, y, t_final, RTOL, ATOL, ITOL, solout, IOUT, WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
	if (IDID == 1) then
		write(*,*) 'Process ended successfully. See RESULT file.'
	else 
		write(*,*) 'Error code:', IDID
	endif
end
