SUBDIRS = src

if ENABLE_PYTHON
SUBDIRS += pytools
endif

ACLOCAL_AMFLAGS = -I m4

if ENABLE_PYTHON
# make verify_cluster_algorithm.py an executable script
install-exec-hook:
	chmod 755 $(bindir)/pytools/verify_cluster_algorithm.py
	chmod 755 $(bindir)/pytools/plot_scaling.py
endif
