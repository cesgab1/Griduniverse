    module GBoost
    ! Early-universe boost of Newton's constant: G(a)/G0 = 1 + delta/(1 + (a/a_t)^4)
    ! (grid-stretching phase of the River-Compactness framework). delta = 0 gives standard CAMB.
    use precision
    implicit none
    real(dl) :: gboost_delta = 0._dl, gboost_at = 1.e-4_dl
    contains
    function Gfac(a)
    real(dl), intent(in) :: a
    real(dl) :: Gfac
    if (gboost_delta == 0._dl) then
        Gfac = 1._dl
    else
        Gfac = 1._dl + gboost_delta/(1._dl + (a/gboost_at)**4)
    end if
    end function Gfac
    function dlnGfac_dlna(a)
    real(dl), intent(in) :: a
    real(dl) :: dlnGfac_dlna, x
    if (gboost_delta == 0._dl) then
        dlnGfac_dlna = 0._dl
    else
        x = (a/gboost_at)**4
        dlnGfac_dlna = -4._dl*gboost_delta*x/(1._dl + x)**2/Gfac(a)
    end if
    end function dlnGfac_dlna
    subroutine set_gboost(delta, at) bind(C, name='set_gboost')
    use iso_c_binding
    real(c_double), value :: delta, at
    gboost_delta = delta
    gboost_at = at
    end subroutine set_gboost
    end module GBoost
