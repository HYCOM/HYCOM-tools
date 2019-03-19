#
# Makefile for HYCOM transport sections processing
#
# --- usage: make transport ARCH=intelGF >& Make_transport.log
#
# --- Tunable parameters in ../../config/$(ARCH)_setup
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
	@exit 1
.f:
	@echo "Must have an explicit rule for" $*
	@exit 1
.c:
	@echo "Must have an explicit rule for" $*
	@exit 1

include ../../config/$(ARCH)_setup

#
# ---------------------------------------------------------------------
# Standard part of the Makefile starts here
# ---------------------------------------------------------------------
#

all:
	/bin/csh Make_all.csh

clean:
	/bin/rm *.o *.a *.mod M*log *tspt *port *port2 *port3 *_lm *vel *mn *mn2 *mn3 *2p0

MODS =	mod_trans.o mod_xc.o mod_za.o wtime.o
OBJS =	blkin.o geopar.o getdat.o zh.o

#
# --- executables.
#

baro_vel:        $(MODS) baro_vel.o              $(OBJS)
	$(LD) $(LDFLAGS) baro_vel.o      $(MODS) $(OBJS) -o baro_vel

baro_vel_mn:     $(MODS) baro_vel_mn.o           $(OBJS)
	$(LD) $(LDFLAGS) baro_vel_mn.o   $(MODS) $(OBJS) -o baro_vel_mn

transport:       $(MODS) transport.o             $(OBJS)
	$(LD) $(LDFLAGS) transport.o     $(MODS) $(OBJS) -o transport

transport2:      $(MODS) transport2.o            $(OBJS)
	$(LD) $(LDFLAGS) transport2.o    $(MODS) $(OBJS) -o transport2

transport3:      $(MODS) transport3.o            $(OBJS)
	$(LD) $(LDFLAGS) transport3.o    $(MODS) $(OBJS) -o transport3

transport3_lm:   $(MODS) transport3_lm.o         $(OBJS)
	$(LD) $(LDFLAGS) transport3_lm.o $(MODS) $(OBJS) -o transport3_lm

transp_mn:       $(MODS) transp_mn.o             $(OBJS)
	$(LD) $(LDFLAGS) transp_mn.o     $(MODS) $(OBJS) -o transp_mn

transp_mn2:      $(MODS) transp_mn2.o            $(OBJS)
	$(LD) $(LDFLAGS) transp_mn2.o    $(MODS) $(OBJS) -o transp_mn2

transp_mn3:      $(MODS) transp_mn3.o            $(OBJS)
	$(LD) $(LDFLAGS) transp_mn3.o    $(MODS) $(OBJS) -o transp_mn3

transp_mn_2p0:   $(MODS) transp_mn_2p0.o         $(OBJS)
	$(LD) $(LDFLAGS) transp_mn_2p0.o $(MODS) $(OBJS) -o transp_mn_2p0

meantspt:                meantspt.o              blkin.o zh.o
	$(LD) $(LDFLAGS) meantspt.o              blkin.o zh.o -o meantspt

mergetspt:               mergetspt.o             blkin.o zh.o
	$(LD) $(LDFLAGS) mergetspt.o             blkin.o zh.o -o mergetspt

#
# --- explicit dependencies for each .[fF] file.
#

baro_vel.o:      baro_vel.f      $(MODS)
baro_vel_mn.o:   baro_vel_mn.f   $(MODS)
transport.o:     transport.f  $(MODS)
transport2.o:    transport2.f $(MODS)
transport3.o:    transport3.f $(MODS)
transp_mn.o:     transp_mn.f  $(MODS)
transp_mn2.o:    transp_mn2.f $(MODS)
transp_mn3.o:    transp_mn3.f $(MODS)
transp_mn_2p0.o: transp_mn_2p0.f $(MODS)
geopar.o:        geopar.f     $(MODS)
getdat.o:        getdat.f     $(MODS)

blkin.o:     blkin.f
meantspt.o:  meantspt.f
mergetspt.o: mergetspt.f
wtime.o:     wtime.F

mod_trans.o: mod_trans.F mod_xc.o mod_za.o
mod_xc.o:    mod_xc.F   
mod_za.o:    mod_za.F    mod_xc.o
mod_zb.o:    mod_zb.F

zh.o:        zh.F zh_sun.f zh_t3e.f zh_ymp.f
