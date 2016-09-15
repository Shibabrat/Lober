bin_SCRIPTS = lober-auto
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = lober-auto.sh

lober-auto: $(srcdir)/lober-auto.sh
	rm -f lober-auto
	echo "#! " $(SH) > lober-auto
	cat $(srcdir)/lober-auto.sh >> lober-auto
	chmod ugo+x lober-auto
