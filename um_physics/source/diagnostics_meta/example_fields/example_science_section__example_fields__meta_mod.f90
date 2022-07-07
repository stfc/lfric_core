!-------------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!
!               Sample Module Diagnostic Meta File
!
! This is a sample file that illustrates how you define meta data for a field
! THIS IS TO BE UPDATED WITH DETAILS ABOUT WHERE THESE FILES ARE TO BE LOCATED,
! ONCE THIS IS DECIDED.
!
! (*/source/diagnostics_meta/*_meta_mod.f90)
!
!-------------------------------------------------------------------------------
module example_science_section__example_fields__meta_mod

  use diagnostics_mod,                only: field_meta_data_type
  use constants_mod,                  only: real_type, r_def, i_def, &
                                            str_short, str_def
  !> Only import the dimensions that you will actually be using
  use vertical_dimensions_mod,        only: model_height_dimension, &
                                            model_depth_dimension, &
                                            fixed_height_dimension
  use non_spatial_dimension_mod,      only: non_spatial_dimension_type, &
                                            NUMERICAL, &
                                            CATEGORICAL
  use misc_meta_data_mod,             only: misc_meta_data_type
  use field_synonym_mod,              only: field_synonym_type
  !> Only import the function spaces that you will actually be using
  use fs_continuity_mod,              only: W2H, W3, Wtheta
  !> Only import the time steps that you will actually be using
  use time_step_enum_mod,             only: STANDARD_TIMESTEP
  !> Only import the interpolation methods that you will actually be using
  use interpolation_enum_mod,         only: BILINEAR
  !> Only import the levels that you will actually be using
  use levels_enum_mod,                only: TOP_WET_LEVEL, &
                                            BOTTOM_SOIL_LEVEL, &
                                            TOP_SOIL_LEVEL, &
                                            BOTTOM_ATMOSPHERIC_LEVEL, &
                                            TOP_ATMOSPHERIC_LEVEL
  use positive_enum_mod,              only: POSITIVE_UP, POSITIVE_DOWN
  use field_synonyms_enum_mod,        only: AMIP, GRIB, CF, CMIP6, STASH

  implicit none

  private

  type, public :: example_science_section__example_fields__meta_type

    !> Declare the name of your fields here
    type(field_meta_data_type), public :: &
      eastward_wind, &
      rate_of_increase_rain_mass_due_to_autoconv_from_liquid_cloud, &
      air_potential_temperature, &
      moisture_content_of_soil_layer, &
      surface_altitude, &
      air_temperature_over_tiles, &
      low_type_cloud_area_fraction
      character(str_def) :: name = "example_science_section__example_fields"

    end type example_science_section__example_fields__meta_type

  interface example_science_section__example_fields__meta_type
    module procedure example_science_section__example_fields__meta_constructor
  end interface

contains

  !>@brief Creates field_meta_data_type objects for a specific section of science
  function example_science_section__example_fields__meta_constructor() result(self)
    implicit none

    type(example_science_section__example_fields__meta_type) :: self

    !> Example field using a height vertical dimension using model levels
    !> If no arguments are present, it will default to top atmospheric level
    !> and bottom atmospheric level. This field also uses a standard name,
    !> which is optional
    self%eastward_wind = field_meta_data_type(&
      unique_id = "example_fields__eastward_wind", &
      units = "m s-1", &
      function_space = W2H, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "u component of wind on u pts on native c grid.", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension( &
              bottom = BOTTOM_ATMOSPHERIC_LEVEL, &
              top = TOP_ATMOSPHERIC_LEVEL), &
      standard_name = "eastward_wind", &
      synonyms = [ &
            field_synonym_type(STASH, "2"),& !> literally this stash code or approx - let the user know
            field_synonym_type(AMIP, "ua"),&
            field_synonym_type(GRIB, "33 E131"),&
            field_synonym_type(CF, "eastward_wind"),&
            field_synonym_type(CMIP6, "ua")&
        ],&
      misc_meta_data = [misc_meta_data_type("positive","eastwards")])


    !> Example of a field using a model height dimension and supplying one
    !> of its arguments
    self%rate_of_increase_rain_mass_due_to_autoconv_from_liquid_cloud = field_meta_data_type(&
      unique_id = "example_fields__rate_of_increase_of_rain_mass_due_to_"// &
                  "autoconversion_from_liquid_cloud", &
      long_name = "rate_of_increase_rain_mass_due_to_autoconv_from_liquid_cloud", &
      units = "kg kg-1 s-1", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "This is a microphysical process transfer rate. "// &
                    "Outputting all microphysics", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension( &
              bottom = BOTTOM_ATMOSPHERIC_LEVEL, &
              top = TOP_WET_LEVEL))

    !> This field uses a different function space
    self%air_potential_temperature = field_meta_data_type(&
      unique_id = "example_fields__air_potential_temperature", &
      units = "K", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "Potential temperature on p points on native c grid.", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_height_dimension( &
              bottom = BOTTOM_ATMOSPHERIC_LEVEL, &
              top = TOP_ATMOSPHERIC_LEVEL), &
      standard_name = "air_potential_temperature",&
      synonyms = [ &
            field_synonym_type(AMIP, "theta"),&
            field_synonym_type(CF, "air_potential_temperature"),&
            field_synonym_type(GRIB, "13")&
            ])

    !> Example field using a depth vertical dimension using model levels
    !> If no arguments are present, it will default to top soil level and
    !> bottom soil level
    self%moisture_content_of_soil_layer = field_meta_data_type(&
      unique_id = "example_fields__moisture_content_of_soil_layer", &
      units = "kg m-2", &
      function_space = W3, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "Total (frozen+unfrozen) soil moisture content in a "// &
                    "soil layer (kg/m2). This example field will be dealt"// &
                    " with using negative height values (Phase 3)", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = model_depth_dimension( &
              bottom = BOTTOM_SOIL_LEVEL, &
              top = TOP_SOIL_LEVEL), &
      standard_name = "mass_content_of_water_in_soil_layer",&
      synonyms = [ &
            field_synonym_type(CF, "mass_content_of_water_in_soil_layer"),&
            field_synonym_type(STASH, "9")&
            ])

    !> Example field using a height vertical dimension on fixed levels
    !> The level definition should be passed as an array of
    !> floating point numbers
    self%surface_altitude = field_meta_data_type(&
      unique_id = "example_fields__surface_altitude", &
      units = "m", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = 'The surface called "surface" means the lower boundary '// &
                    'of the atmosphere. Altitude is the (geometric) height '// &
                    'above the geoid, which is the reference geopotential '// &
                    'surface. The geoid is similar to mean sea level.', &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension( &
                                      level_definition = REAL([0.0], r_def)), &
      non_spatial_dimension = [non_spatial_dimension_type( &
              dimension_name = "test_axis_non_spatial_dimension", &
              dimension_category = NUMERICAL, &
              help_text = "Test axis non-spatial dimension help text", &
              non_spatial_units = "test_unit_1"), &
      non_spatial_dimension_type( &
              dimension_name = "test_axis_non_spatial_dimension_2", &
              dimension_category = CATEGORICAL, &
              help_text = "Test axis non-spatial dimension help text_2", &
              label_definition = [character(str_short) :: 'A','B','C','D','E'], &
              non_spatial_units = "test_unit_2")], &
      standard_name = "surface_altitude", &
      synonyms = [ &
            field_synonym_type(AMIP, "orog"),&
            field_synonym_type(CF, "surface_altitude"),&
            field_synonym_type(CMIP6, "orog"),&
            field_synonym_type(STASH, "33")&
            ])

    !> Another example of a field using a height vertical
    !> dimension on fixed levels.
    !> This field will require an additional non-spatial dimension to
    !> define the tiles. This is yet to be implemented
    self%air_temperature_over_tiles = field_meta_data_type(&
      unique_id = "example_fields__air_temperature_over_tiles", &
      units = "K", &
      function_space = Wtheta, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "1.5M TEMPERATURE OVER TILES", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension( &
              level_definition = real([1.5], r_def)), &
      non_spatial_dimension = [non_spatial_dimension_type( &
              dimension_name = "Tiles", &
              dimension_category = CATEGORICAL, &
              help_text = "Tiles help text", &
              label_definition = [character(str_short) :: 'Broadleaf Tree', &
                                                          'Needle Leaf Tree', &
                                                          'C3 Grass', &
                                                          'C4 Grass', &
                                                          'Shrub', &
                                                          'Urban', &
                                                          'Water', &
                                                          'Soil', &
                                                          'Ice']), &
      non_spatial_dimension_type( &
              dimension_name = "mutable_dimension", &
              dimension_category = CATEGORICAL, &
              help_text = "test_mutable_text", &
              non_spatial_units = "test_mutable_unit"), &
      non_spatial_dimension_type( &
              dimension_name = "mutable_dimension_no_units", &
              dimension_category = CATEGORICAL, &
              help_text = "test_mutable_text")], &
      standard_name = "air_temperature", &
      synonyms = [ &
            field_synonym_type(STASH, "3328")&
            ])

    !> Example of a field using a fixed height dimension
    self%low_type_cloud_area_fraction = field_meta_data_type(&
      unique_id = "example_fields__low_type_cloud_area_fraction", &
      units = "1", &
      function_space = W3, &
      order = 0, &
      io_driver = "", &
      trigger = "__checksum: true;", &
      description = "Low cloud amount. Phase 3", &
      data_type = REAL_TYPE, &
      time_step = STANDARD_TIMESTEP, &
      recommended_interpolation = BILINEAR, &
      packing = 0, &
      vertical_dimension = fixed_height_dimension( &
              level_definition = REAL([111.0, 1949.0], r_def)), &
      standard_name = "low_type_cloud_area_fraction",&
      synonyms = [ &
            field_synonym_type(STASH, "9203")&
            ])

  end function example_science_section__example_fields__meta_constructor
end module example_science_section__example_fields__meta_mod