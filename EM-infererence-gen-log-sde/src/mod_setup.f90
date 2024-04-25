! Created by saul on 11/17/23.
module mod_setup
    

contains
    !subroutine load_parameters(filename)
    subroutine load_parameters()
        character (100)            :: filename
        
        filename = "./data/config.json"
        print *, "The namefile is: ", trim(filename)        
        !
        !if (json%failed()) stop
        !
    end subroutine load_parameters
end module  mod_setup