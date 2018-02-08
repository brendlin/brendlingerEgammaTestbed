################################################################################
# Package: brendlingerEgammaTestbed
################################################################################

# Declare the package name:
atlas_subdir( brendlingerEgammaTestbed )

# Declare the package dependencies:
atlas_depends_on_subdirs( PUBLIC
                          HGamCore/HGamAnalysisFramework
                          )

# Find the needed external(s):
find_package( ROOT COMPONENTS Core RIO Hist Tree )

# build a CINT dictionary for the library
atlas_add_root_dictionary ( brendlingerEgammaTestbedLib brendlingerEgammaTestbedCintDict
                            ROOT_HEADERS brendlingerEgammaTestbed/RadiativeZSelection.h Root/LinkDef.h
                            EXTERNAL_PACKAGES ROOT
                            )

# build a shared library
atlas_add_library( brendlingerEgammaTestbedLib
                   brendlingerEgammaTestbed/*.h Root/*.cxx ${brendlingerEgammaTestbedCintDict}
                   PUBLIC_HEADERS brendlingerEgammaTestbed
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
                   HGamAnalysisFrameworkLib
                   )

# Install files from the package:
atlas_install_data( data/* )

# Install scripts from the package:
# atlas_install_scripts( scripts/*.sh scripts/*.py )

atlas_add_executable( runRadiativeZSelection util/runRadiativeZSelection.cxx
                      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                      LINK_LIBRARIES ${ROOT_LIBRARIES} brendlingerEgammaTestbedLib )
