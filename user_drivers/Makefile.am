# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= thin_film_brinkman1D_refineable \
thin_film_brinkman2D_triangle \
thin_film_brinkman2D_ellipse_structured
#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
thin_film_brinkman1D_refineable_SOURCES = thin_film_brinkman1D_refineable.cc
thin_film_brinkman2D_triangle_SOURCES = thin_film_brinkman2D_triangle.cc
thin_film_brinkman2D_ellipse_structured_SOURCES = thin_film_brinkman2D_ellipse_structured.cc elliptic_integral.cc

# Required libraries:
thin_film_brinkman1D_refineable_LDADD = -L@libdir@ -lthin_film_brinkman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
thin_film_brinkman2D_triangle_LDADD = -L@libdir@ -lthin_film_brinkman -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
thin_film_brinkman2D_ellipse_structured_LDADD = -L@libdir@ -lthin_film_brinkman -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS += -I@includedir@
