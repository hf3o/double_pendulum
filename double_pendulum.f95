program double_pendulum
  implicit none
  real(8) :: t, dt, t_end
  real(8) :: theta1, theta2, theta1_dot, theta2_dot
  real(8) :: m1, m2, l1, l2, g
  real(8) :: x2, y2
  integer :: i
  open(unit=10, file="trajectory_data.dat", status="replace")

  ! Parametry systemu
  m1 = 1.0
  m2 = 1.0
  l1 = 1.0
  l2 = 1.0
  g = 9.81
  dt = 0.001
  t_end = 10.0

  ! Początkowe warunki
  theta1 = 0.1
  theta2 = 0.15
  theta1_dot = 0.0
  theta2_dot = 0.0

  ! Nagłówki w pliku wyjściowym
  write(10, '(A)') "# x2 y2"

  do i = 1, int(t_end / dt)
     t = (i-1) * dt

     ! Oblicz współrzędne końcówki wahadła
     x2 = l1 * sin(theta1) + l2 * sin(theta2)
     y2 = -l1 * cos(theta1) - l2 * cos(theta2)

     ! Sprawdź, czy wartości są NaN
     if (isnan(x2) .or. isnan(y2)) then
         write(*,*) "NaN detected at t =", t
         stop
     end if

     ! Zapisz dane do pliku, oddzielone spacjami
     write(10, '(F10.6,1X,F10.6)') x2, y2

     ! Aktualizacja kątów
     call update_theta(theta1, theta2, theta1_dot, theta2_dot, m1, m2, l1, l2, g, dt)
  end do

  close(10)
contains

  subroutine update_theta(theta1, theta2, theta1_dot, theta2_dot, m1, m2, l1, l2, g, dt)
    real(8), intent(inout) :: theta1, theta2, theta1_dot, theta2_dot
    real(8), intent(in) :: m1, m2, l1, l2, g, dt
    real(8) :: delta_theta, theta1_ddot, theta2_ddot

    delta_theta = theta2 - theta1
    theta1_ddot = (m2 * l2 * theta2_dot**2 * sin(delta_theta) + m2 * g * sin(theta2) * cos(delta_theta) + m2 * l2 * theta1_dot**2 * sin(delta_theta) * cos(delta_theta) - (m1 + m2) * g * sin(theta1)) / ((m1 + m2) * l1 - m2 * l1 * cos(delta_theta)**2)
    theta2_ddot = (-m2 * l2 * theta2_dot**2 * sin(delta_theta) * cos(delta_theta) + (m1 + m2) * g * sin(theta1) * cos(delta_theta) - (m1 + m2) * l1 * theta1_dot**2 * sin(delta_theta) - (m1 + m2) * g * sin(theta2)) / (l2 / l1 * (m1 + m2) * sin(delta_theta))

    theta1_dot = theta1_dot + theta1_ddot * dt
    theta2_dot = theta2_dot + theta2_ddot * dt
    theta1 = theta1 + theta1_dot * dt
    theta2 = theta2 + theta2_dot * dt
  end subroutine

  logical function isnan(x)
    real(8), intent(in) :: x
    isnan = (x /= x) ! Wartość NaN nie jest równa samej sobie
  end function isnan

end program double_pendulum

