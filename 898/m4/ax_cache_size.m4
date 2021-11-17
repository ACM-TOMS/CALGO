# ===========================================================================
#             http://autoconf-archive.cryp.to/ax_cache_size.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CACHE_SIZE
#
# DESCRIPTION
#
#   Find L1 and L2 caches size by reading the corresponding file on UNIX or
#   by requesting cpuid. The results are available in the defines
#   CPU_L1_CACHE and CPU_L2_CACHE.
#
#   This macro depends on AX_GCC_X86_CPUID, AC_PROG_SED, AX_COUNT_CPUS, and
#   AX_CPU_VENDOR.
#
# LAST MODIFICATION
#
#   2008-04-12
#
# COPYLEFT
#
#   Copyright (c) 2008 Christophe Tournayre <turn3r@users.sourceforge.net>
#
# Patched by:
#
#   Copyright (c) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#   Copyright (c) 2008 Arnaud Bergeron <abergeron@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_CACHE_SIZE],
[
  AC_REQUIRE([AC_PROG_SED])
  AC_REQUIRE([AX_COUNT_CPUS])
  AC_REQUIRE([AX_GCC_X86_CPUID])
  AC_REQUIRE([AX_CPU_VENDOR])

  AX_COUNT_CPUS
  AX_CPU_VENDOR

    ax_l1_size=unknown
    ax_l2_size=unknown

    #Check if the variable is present
    if test -e /sys/devices/system/cpu/cpu0/cache/index0/size; then
      for ncpu in `seq 0 $(($CPU_COUNT-1))`; do
        for idx in `seq 0 3`; do
          if test -e /sys/devices/system/cpu/cpu$ncpu/cache/index$idx/size ; then
            level=`cat /sys/devices/system/cpu/cpu$ncpu/cache/index$idx/level`
            size=`cat /sys/devices/system/cpu/cpu$ncpu/cache/index$idx/size`
            eval CPU$ncpu\_L$level\_CACHE="$size"
          fi
        done
      done

      # This part can (must !!!) be optimized, because we know all caches per proc but
      # we only take care about the first proc
      ax_l1_size=$CPU0_L1_CACHE
      ax_l2_size=$CPU0_L2_CACHE

    else
      if test "x$ax_cpu_vendor" != "xUnknown"; then
        #Or use CPUID
        if test "x$ax_cpu_vendor" != "xIntel"; then
          AX_GCC_X86_CPUID(0x80000005) # For L1 cache (not available on intel !!!)
        
          l1_hexval=$(( 16#`echo $ax_cv_gcc_x86_cpuid_0x80000005 | cut -d ":" -f 4`))
          ax_l1_size=$(($l1_hexval >> 24))
        fi
        
        AX_GCC_X86_CPUID(0x80000006) # For L2 cache
        
        l2_hexval=$(( 16#`echo $ax_cv_gcc_x86_cpuid_0x80000006 | cut -d ":" -f 3`))
        ax_l2_size=$(($l2_hexval >> 16))
      else
        #Or use sysctl
        if test "`/usr/sbin/sysctl -n hw.l2cachesize 2>/dev/null`"; then
          ax_l2_size=$((`/usr/sbin/sysctl -n hw.l2cachesize` / 1024))
        fi
        
        if test "`/usr/sbin/sysctl -n hw.l1dcachesize 2>/dev/null`"; then
          ax_l1_size=$((`/usr/sbin/sysctl -n hw.l1dcachesize` / 1024))
        fi
        if test "`/usr/sbin/sysctl -n hw.l1cachesize 2>/dev/null`"; then
          ax_l1_size=$((`/usr/sbin/sysctl -n hw.l1cachesize` / 1024))
        fi
      fi
    fi

    # Keep only digits if there is a unit (ie 1024K -> 1024) and convert in Bytes
    AC_MSG_CHECKING(the L1 cache size)
    ax_l1_size=`echo $ax_l1_size | $SED 's/\([[0-9]]\)[[A-Za-z]]$/\1/g'`
    ax_l1_size=$(($ax_l1_size*1024))
    AC_MSG_RESULT( $ax_l1_size Bytes)

    AC_MSG_CHECKING(the L2 cache size)
    ax_l2_size=`echo $ax_l2_size | $SED 's/\([[0-9]]\)[[A-Za-z]]$/\1/g'`
    ax_l2_size=$(($ax_l2_size*1024))
    AC_MSG_RESULT( $ax_l2_size Bytes)

    AC_DEFINE_UNQUOTED([CPU_L1_CACHE], ${ax_l1_size}, [L1 cache size (in Bytes)])
    AC_DEFINE_UNQUOTED([CPU_L2_CACHE], ${ax_l2_size}, [L2 cache size (in Bytes)])
])
