#AC_SEARCH_LHAPDF(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_HOPPET],[
  AC_ARG_WITH([hoppet], AC_HELP_STRING(--with-hoppet, [path to HOPPET library and header files]))

  ## Use a specified --with-lhapdf arg to set basic paths, if provided
  HOPPETCONFIG_PATH=$PATH
  if test -e "$with_lhapdf"; then
    HOPPETCONFIG_PATH=$with_lhapdf/bin:$HOPPETCONFIG_PATH
    HOPPETPATH="$with_lhapdf"
    HOPPETINCPATH="$HOPPETPATH/include"
    HOPPET_CPPFLAGS="-I$HOPPETINCPATH"
    HOPPET_CXXFLAGS=""
  fi

  ## Try to do better, using the lhapdf-config script
  AC_PATH_PROG(HOPPETCONFIG, hoppet-config, [], [$HOPPETCONFIG_PATH])
  if test -x "$HOPPETCONFIG"; then
    AC_MSG_NOTICE(Using $HOPPETCONFIG to find HOPPET flags)
    HOPPETPATH=`$HOPPETCONFIG --prefix`
    HOPPETINCPATH="$HOPPETPATH/include"
    HOPPET_CPPFLAGS=`$HOPPETCONFIG --cxxflags`
    HOPPET_CXXFLAGS=`$HOPPETCONFIG --cxxflags`
    HOPPET_LDFLAGS=`$HOPPETCONFIG --libs`
    HOPPET_FFFLAGS=`$HOPPETCONFIG --fflags`
  fi

  ## If it's worked, propagate the variables and execute success arg
  if test -e "$HOPPETPATH"; then
    ## Otherwise  execute the fail arg
    AC_SUBST([HOPPETPATH])
    AC_SUBST([HOPPETINCPATH])
    AC_SUBST([HOPPET_CPPFLAGS])
    AC_SUBST([HOPPET_CXXFLAGS])
    AC_SUBST([HOPPET_LDFLAGS])
    AC_SUBST([HOPPET_FFFLAGS])
    AM_CONDITIONAL([WITH_HOPPET], true)
    AM_CONDITIONAL([WITH_HOPPETLIB], true)
    AM_CONDITIONAL([WITH_HOPPETINC], true)
    AM_CONDITIONAL([WITHOUT_HOPPET], false)
    AM_CONDITIONAL([WITHOUT_HOPPETLIB], false)
    AM_CONDITIONAL([WITHOUT_HOPPETINC], false)
    AC_MSG_NOTICE([HOPPET include path is $HOPPETINCPATH])
    AC_MSG_NOTICE([HOPPET CPPFLAGS is $HOPPET_CPPFLAGS])
    AC_MSG_NOTICE([HOPPET CXXFLAGS is $HOPPET_CXXFLAGS])
    AC_MSG_NOTICE([HOPPET LDFLAGS is $HOPPET_LDFLAGS])
    $1
  else
    AM_CONDITIONAL([WITH_HOPPET], false)
    AM_CONDITIONAL([WITH_HOPPETLIB], false)
    AM_CONDITIONAL([WITH_HOPPETINC], false)
    AM_CONDITIONAL([WITHOUT_HOPPET], true)
    AM_CONDITIONAL([WITHOUT_HOPPETLIB], true)
    AM_CONDITIONAL([WITHOUT_HOPPETINC], true)
    $2
  fi
])
