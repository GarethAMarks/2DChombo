# -*- Mode: Makefile -*-


### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.
###$(GRCHOMBO_SOURCE)/../../GRChombo/Source/ApparentHorizonFinder

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

DIM = 2

GRCHOMBO_SOURCE = ../../Source
GSL_ROOT = /usr/local/software/spack/2024.2/opt/linux-ubuntu22.04-x86_64_v4/gcc-11.4.0/gsl-2.7.1-bw4k4zxuc5yay66l5mic5saul6geudxr
ebase = Main_ScalarField2D



# Figure out whether to build with Two Punctures
ifdef TWOPUNCTURES_SOURCE

#CXXFLAGS += -I$(GSL_ROOT)/include
#LDFLAGS  += -L$(GSL_ROOT)/lib
#LDLIBS   += -lgsl -lgslcblas -lm
#INCLUDE_PATH += $(GSL_ROOT)/include

ifneq (,$(wildcard $(TWOPUNCTURES_SOURCE)/TwoPunctures.hpp))
  $(info ================== Using TwoPunctures ==================)
  USE_TWOPUNCTURES := TRUE
  XTRACPPFLAGS += -I$(GSL_ROOT)/include
  XTRALIBFLAGS  = -L$(GSL_ROOT)/lib -lgsl -lgslcblas -lm
  # Note this code requires linking with the GNU Scientific Library
  # In addition to the line below, the include paths to the relevant header
  # files might also need to be passed to the preprocessor/compiler
  #XTRALIBFLAGS = -lgsl
  cxxcppflags = -DUSE_TWOPUNCTURES
else
  $(error TWOPUNCTURES_SOURCE not set correctly)
endif
endif


LibNames = AMRTimeDependent AMRTools BoxTools


src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Matter  \
            $(GRCHOMBO_SOURCE)/Cartoon  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/InitialConditions/ScalarFields \
            $(GRCHOMBO_SOURCE)/InitialConditions/BosonStars  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles  \
            $(GRCHOMBO_SOURCE)/InitialConditions/ThinShell  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BosonStars/SingleStar \
            $(GRCHOMBO_SOURCE)/ApparentHorizonFinder  \
            $(GRCHOMBO_SOURCE)/BlackHoles \

ifeq ($(USE_TWOPUNCTURES),TRUE)
  src_dirs += $(TWOPUNCTURES_SOURCE)
endif


include $(CHOMBO_HOME)/mk/Make.test



