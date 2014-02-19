!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! Collect together the various parts of LFRic infrastructure into a
! high-level meta-module.
!-------------------------------------------------------------------------------

module lfric
  use constants_mod
  use function_space_mod, only: function_space_type
  use field_mod,          only: field_type
  use kernel_mod,         only: kernel_type
  use psy_mod,            only: psy_type
end module lfric
