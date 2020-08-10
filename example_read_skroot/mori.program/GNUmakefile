#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
#SKOFL_ROOT = /skofl
#
#Note: it should change into csh so that the makefile can compile

include    $(SKOFL_ROOT)/config.gmk

LOCAL_INC	= 
LOCAL_LIBS	= -lsklowe_7.0 -lsollib_4.0 -lwtlib_5.1 -lbonsai_3.3 -lstmu -lska
#NOTE: stmu & ska are ATMPD lib, see /home/atmpd/skrep/09e

# the way this is working is that within config.gmk there is:
#LDLIBS =  $(LDOPTS) $(OBJFJT) $(LOCAL_LIBS) $(LIBS) $(LIBS) $(CLIB) $(ROOTLIBS) \
#	   $(SITE_LIBRARIES) $(FORTRANLIBRARIES) $(SITE_LDLIBS)
# i.e., ${LDLIBS}, which is basically a compiled-in location for dynamically linked libraries,
# has been pre-set to include the list ${LOCAL_LIBS}, which we're then providing.
#
# similarly, it defines:
#SITE_INCLUDES  = $(LOCAL_INC) -I$(CERN_ROOT)/include -I$(CERN_ROOT)/include/cfortran $(ROOTCFLAGS)
# which in turn gets put into: ${CINCLUDES} and ${CXXINCLUDES}
# so that ${LOCAL_INC} lets us append to this list.

# see https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
# for a list of gnu make implicit vars (or  http://web.mit.edu/gnu/doc/html/make_10.html#SEC89)

#
#  Objects
#

OBJS   = dump_rootfile

all: $(OBJS) 
	$(RM) *.o *~

clean: 
	$(RM) *.o *~ core fort.* $(OBJS)

install.exec: 

