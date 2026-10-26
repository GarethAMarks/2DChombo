# -*- Mode: Makefile -*-
# Single-box unit test of the 2D Dirac code (see box_test_dirac.cpp). Set CHOMBO_HOME as for the examples.

DIM = 2

GRCHOMBO_SOURCE = ../../../../Source

ebase = box_test_dirac

LibNames = AMRTimeDependent AMRTools BoxTools

# the example's UserVariables.hpp / DiagnosticVariables.hpp
XTRACPPFLAGS += -I../..

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/Matter  \
            $(GRCHOMBO_SOURCE)/Cartoon  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator  \
            $(GRCHOMBO_SOURCE)/InitialConditions/DiracStars \
            $(GRCHOMBO_SOURCE)/InitialConditions/BosonStars
include $(CHOMBO_HOME)/mk/Make.test
