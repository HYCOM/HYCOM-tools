#
set echo
#
foreach f ( README* *.csh Makefile *.h *.f *.F )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f ~/HYCOM-tools/MPI/meanstd/src
end
