SUBDIRS = src cpp pytools

ACLOCAL_AMFLAGS = -I m4

# make verify_cluster_algorithm.py an executable script
install-exec-hook:
	chmod 755 $(bindir)/pytools/verify_cluster_algorithm.py
	chmod 755 $(bindir)/pytools/plot_scaling.py
