#
set echo
#
#foreach f ( baro_vel.f baro_vel_mn.f meantspt.f mergetspt.f transp_mn.f transp_mn2.f transp_mn3.f transp_mn_2p0.f transport.f transport2.f transport3.f )
foreach f ( baro_vel.f baro_vel_mn.f meantspt.f mergetspt.f transp_mn_2p0.f transp_mn3.f transport3.f transport3_lm.f )
  mv $f ${f}+02
  sed -e "s/\*80/*240/g" ${f}+02 >! $f
end
