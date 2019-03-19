#!/bin/csh 
#set echo
#
# --- Usage:  ./Make_clean.csh >& Make_clean.log
#
# --- set ARCH to the correct value for this machine.
#
source Make_all.src
make clean ARCH=${ARCH}
