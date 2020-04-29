AC_PREREQ([2.63])
AC_INIT(ScalarEMT,[v1.0],[joseph.lee@ed.ac.uk],[ScalarEMT])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([scalar_EMT.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [try this for a non-standard install prefix of Hadrons])],
    [PATH="$with_hadrons/bin$PATH_SEPARATOR$PATH"]
    [CXXFLAGS="$CXXFLAGS -I$with_hadrons/include"]
    [LDFLAGS="$LDFLAGS -L$with_hadrons/lib"])

AC_CHECK_PROG([HADCONF],[hadrons-config],[yes])
if test x"$HADCONF" != x"yes" ; then
    AC_MSG_ERROR([hadrons-config not found])
fi
CXXFLAGS="$CXXFLAGS `hadrons-config --cxxflags`"
LDFLAGS="$LDFLAGS `hadrons-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS=" -lHadrons $LIBS `hadrons-config --libs`"

AC_PROG_CXX
AC_LANG([C++])

AC_MSG_CHECKING([that a minimal Hadrons program compiles]);
AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
    #include <Hadrons/Application.hpp>

    using namespace Grid;
    using namespace Hadrons;

    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Grid_finalize();

        Application application;

        application.run();

        return 0;
    }

    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Failed to compile a minimal Hadrons program])])

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

