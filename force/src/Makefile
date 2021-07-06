#
# Makefile for HYCOM wind interpolation program
#
# --- usage: make zero ARCH=sun >& Make_zero
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
	/bin/rm -f *.o *.a *.mod M*log M*_?? ?? M*_[k]pc *fc *const *curl *diff *extend *fit *magstress *mm5 *nc *offset *scale *rivers *riv_mon *riv_hf *runoff *corr *zero *monthly *interp *shift *add *climo *diurnal *table *sal *core *barom *margin *wd stoch

MODS = mod_xc.o mod_za.o wtime.o
LIBS = interp.o zh.o
MODP = mod_xc.o mod_za.o wtime.o mod_scrip.o
LIBP = interp.o zh.o scrip.o

# select one of the programs to make,
#  e.g. make zero ARCH=sun >& Make_zero

ap:                      $(MODS) ap.o           $(LIBS)
	$(LD) $(LDFLAGS) ap.o           $(MODS) $(LIBS) -o ap

time_interp:             $(MODS) time_interp.o  $(LIBS)
	$(LD) $(LDFLAGS) time_interp.o  $(MODS) $(LIBS) -o time_interp

time_shift:              $(MODS) time_shift.o  $(LIBS)
	$(LD) $(LDFLAGS) time_shift.o  $(MODS) $(LIBS) -o time_shift

aphf_add:                $(MODS) aphf_add.o     $(LIBS)
	$(LD) $(LDFLAGS) aphf_add.o     $(MODS) $(LIBS) -o aphf_add

aphf_climo:              $(MODS) aphf_climo.o   $(LIBS)
	$(LD) $(LDFLAGS) aphf_climo.o   $(MODS) $(LIBS) -o aphf_climo

aphf_diurnal:            $(MODS) aphf_diurnal.o $(LIBS)
	$(LD) $(LDFLAGS) aphf_diurnal.o $(MODS) $(LIBS) -o aphf_diurnal

aphf_extend:             $(MODS) aphf_extend.o  $(LIBS)
	$(LD) $(LDFLAGS) aphf_extend.o  $(MODS) $(LIBS) -o aphf_extend

aphf_meanfit:            $(MODS) aphf_meanfit.o $(LIBS)
	$(LD) $(LDFLAGS) aphf_meanfit.o $(MODS) $(LIBS) -o aphf_meanfit

aphf_monthly:            $(MODS) aphf_monthly.o $(LIBS)
	$(LD) $(LDFLAGS) aphf_monthly.o $(MODS) $(LIBS) -o aphf_monthly

aphf_offset:             $(MODS) aphf_offset.o  $(LIBS)
	$(LD) $(LDFLAGS) aphf_offset.o  $(MODS) $(LIBS) -o aphf_offset

aphf_margin:             $(MODS) aphf_margin.o  $(LIBS)
	$(LD) $(LDFLAGS) aphf_margin.o  $(MODS) $(LIBS) -o aphf_margin

aphf_scale:              $(MODS) aphf_scale.o   $(LIBS)
	$(LD) $(LDFLAGS) aphf_scale.o   $(MODS) $(LIBS) -o aphf_scale

aphf_flcorr:             $(MODS) aphf_flcorr.o  $(LIBS)
	$(LD) $(LDFLAGS) aphf_flcorr.o  $(MODS) $(LIBS) -o aphf_flcorr

aphf_tacorr:             $(MODS) aphf_tacorr.o  $(LIBS)
	$(LD) $(LDFLAGS) aphf_tacorr.o  $(MODS) $(LIBS) -o aphf_tacorr

kp:                      $(MODS) kp.o           $(LIBS)
	$(LD) $(LDFLAGS) kp.o           $(MODS) $(LIBS) -o kp

kpc:                     $(MODS) kp_const.o   zh.o   
	$(LD) $(LDFLAGS) kp_const.o     $(MODS) zh.o    -o kp_const

kphfc:                   $(MODS) kphf_const.o   zh.o   
	$(LD) $(LDFLAGS) kphf_const.o   $(MODS) zh.o    -o kphf_const

kphf_table:              $(MODS) kphf_table.o   zh.o   
	$(LD) $(LDFLAGS) kphf_table.o   $(MODS) zh.o    -o kphf_table

kp_invbarom:             $(MODS) kp_invbarom.o  $(LIBS)
	$(LD) $(LDFLAGS) kp_invbarom.o  $(MODS) $(LIBS) -o kp_invbarom

stoch:                   $(MODS) stoch.o        $(LIBS)
	$(LD) $(LDFLAGS) stoch.o        $(MODS) $(LIBS) -o stoch

tp:                      $(MODS) tp.o           $(LIBS)
	$(LD) $(LDFLAGS) tp.o           $(MODS) $(LIBS) -o tp

tp_sal:                  $(MODS) tp_sal.o       $(LIBS)
	$(LD) $(LDFLAGS) tp_sal.o       $(MODS) $(LIBS) -o tp_sal

wc:                      $(MODS) w_const.o    zh.o   
	$(LD) $(LDFLAGS) w_const.o      $(MODS) zh.o    -o w_const

wi:                      $(MODS) wi.o           $(LIBS)
	$(LD) $(LDFLAGS) wi.o           $(MODS) $(LIBS) -o wi

wi_curl:                 $(MODS) wi_curl.o      $(LIBS)
	$(LD) $(LDFLAGS) wi_curl.o      $(MODS) $(LIBS) -o wi_curl

wi_ewdnwd:                $(MODS) wi_ewdnwd.o   $(LIBS)
	$(LD) $(LDFLAGS) wi_ewdnwd.o    $(MODS) $(LIBS) -o wi_ewdnwd

wi_xwdywd:                $(MODS) wi_xwdywd.o   $(LIBS)
	$(LD) $(LDFLAGS) wi_xwdywd.o    $(MODS) $(LIBS) -o wi_xwdywd

wi_magstress:         $(MODS) wi_magstress.o    $(LIBS)
	$(LD) $(LDFLAGS) wi_magstress.o $(MODS) $(LIBS) -o wi_magstress

wi_meanfit:              $(MODS) wi_meanfit.o   $(LIBS)
	$(LD) $(LDFLAGS) wi_meanfit.o   $(MODS) $(LIBS) -o wi_meanfit

diff:                    $(MODS) off_diff.o   $(LIBS)
	$(LD) $(LDFLAGS) off_diff.o     $(MODS) $(LIBS) -o off_diff

zero:                    $(MODS) off_zero.o     zh.o
	$(LD) $(LDFLAGS) off_zero.o     $(MODS) zh.o    -o off_zero

pzero:                   $(MODS) pcip_zero.o    zh.o
	$(LD) $(LDFLAGS) pcip_zero.o    $(MODS) zh.o    -o pcip_zero

riv_mon:                 $(MODS) pcip_riv_mon.o $(LIBS)
	$(LD) $(LDFLAGS) pcip_riv_mon.o $(MODS) $(LIBS) -o pcip_riv_mon

riv_hf:                  $(MODS) pcip_riv_hf.o  $(LIBS)
	$(LD) $(LDFLAGS) pcip_riv_hf.o $(MODS)  $(LIBS) -o pcip_riv_hf

runoff:                  $(MODS) pcip_runoff.o  $(LIBS)
	$(LD) $(LDFLAGS) pcip_runoff.o  $(MODS) $(LIBS) -o pcip_runoff

conv_f_mm5:              $(MODS) conv_f_mm5.o   zh.o
	$(LD) $(LDFLAGS) conv_f_mm5.o  $(MODS)  zh.o    -o conv_f_mm5

conv_w_mm5:              $(MODS) conv_w_mm5.o   zh.o
	$(LD) $(LDFLAGS) conv_w_mm5.o  $(MODS)  zh.o    -o conv_w_mm5

force2nc:                            force2nc.o $(MODS) zh.o   
	$(LD) $(LDFLAGS) -o force2nc force2nc.o $(MODS) zh.o    $(EXTRANCDF)

nrl2nc:                              nrl2nc.o   $(MODS) zh.o   
	$(LD) $(LDFLAGS) -o nrl2nc   nrl2nc.o   $(MODS) zh.o    $(EXTRANCDF)

ap_nc:                       $(MODS) ap_nc.o            $(LIBS)
	$(LD) $(LDFLAGS) -o ap_nc    ap_nc.o    $(MODS) $(LIBS) $(EXTRANCDF)

kp_core:                    $(MODS) kp_core.o           $(LIBS)
	$(LD) $(LDFLAGS) -o kp_core kp_core.o   $(MODS) $(LIBS) $(EXTRANCDF)

kp_nc:                       $(MODP) kp_nc.o            $(LIBP)
	$(LD) $(LDFLAGS) -o kp_nc    kp_nc.o    $(MODP) $(LIBP) $(EXTRANCDF)

tp_nc:                       $(MODS) tp_nc.o            $(LIBS)
	$(LD) $(LDFLAGS) -o tp_nc    tp_nc.o    $(MODS) $(LIBS) $(EXTRANCDF)

wi_core:                     $(MODP) wi_core.o          $(LIBP)
	$(LD) $(LDFLAGS) -o wi_core  wi_core.o  $(MODP) $(LIBP) $(EXTRANCDF)

wi_nc:                       $(MODP) wi_nc.o            $(LIBP)
	$(LD) $(LDFLAGS) -o wi_nc    wi_nc.o    $(MODP) $(LIBP) $(EXTRANCDF)

mod_scrip.o: mod_scrip.f

mod_xc.o:    mod_xc.F

mod_za.o:    mod_za.F mod_xc.o

zh.o:        zh.F zh_sun.f zh_t3e.f zh_ymp.f
