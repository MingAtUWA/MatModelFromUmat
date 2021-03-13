program main
    implicit none
    
    real(8)::stress(6)
    real(8)::ddsdde(6, 6)
    real(8)::statev(13)
    real(8)::dstran(6)
    real(8)::props(16)
    integer::error
    
    call sand_hypoplasticity_integration(stress, &
        ddsdde, statev, dstran, props, error)
    
endprogram