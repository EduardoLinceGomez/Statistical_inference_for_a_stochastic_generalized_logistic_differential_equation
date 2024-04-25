program mod_setup_test
    ! use mod_setup
    use  json_module, rk => json_rk
    
    implicit none
    
    character(100) filename
    real(kind=rk)              :: t0, dt, tf, mu
    real(kind=rk), allocatable :: x0(:)
    type(json_file)            :: json
    logical                    :: is_found

    filename = '../data/config.json'
    ! call load_parameters()
    call json%initialize()
    call json%load_file('../data/config.json'); 
    if (json%failed()) stop

    ! Read in the data.
    json_block: block
        call json%get('t0', t0, is_found);
        if (.not. is_found) exit json_block
        call json%get('dt', dt, is_found);
        if (.not. is_found) exit json_block
        call json%get('tf', tf, is_found);
        if (.not. is_found) exit json_block
        call json%get('mu', mu, is_found);
        if (.not. is_found) exit json_block
        call json%get('x0', x0, is_found);
        if (.not. is_found) exit json_block
    end block json_block

    ! Output values.
    if (is_found) then
        print *, t0, dt, tf, mu
        print *, x0
    end if

    ! Clean up.
    call json%destroy()
end program mod_setup_test


! ifort -I./include/ -o example.bin example.f90 ./lib/libjsonfortran.a