#
# Usage: awk -f rivers_daily_total.awk date=$RSDATE days=61 scale=1.0E-6 river=1 river_anom.txt
#        create a /RIVER/ namelist (see pcip_riv_hf.f) for this river TOTAL
#
#   river_anom.txt contains daily input of the form:
#   YYYY-MM-DD TOTAL MONTHLY ANOMALY
#
#   the last total will persist if we run out of days
#
BEGIN   {
	t = 0.0
        d = 0
        }

/^#/	{
	next
	}

	{
	s = substr($0,1,10)
	if (s == date) {
		d = 1
		printf("&RIVER\n")
		printf("  KRIVER = %2d\n",river)
		printf("  TSCALE = %s \n",scale)
		printf("  TRIVER = \n")
		}
	if (d > 0 && d <= days) {
		t = $2
		printf("%s\n",t)
		d = d + 1
		}
	}

END	{
	if (d <= days) {
		for ( i=d; i <= days; i++) printf("%s\n",t)
		}
	printf("/\n")
	}
