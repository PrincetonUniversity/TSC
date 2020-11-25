ifeq (,$(filter object,$(notdir $(CURDIR))))

# The following section determines the machine and code version
# and prepares a directory for object files and modules
#--------------------------------------------------------------

include ../tsc.mk

else

# The machine-independent parts of the makefile are specified in this section.
# Machine-dependent parts should be put in makelib_f90
# ----------------------------------------------------------------------------

VPATH= $(SRCDIR)/source:$(subst tsc_a,modules,$(SRCDIR))

include $(SRCDIR)/makelib_f90

endif
