# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test ! -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

action pair_sdpd_taitwater_isothermal.cpp atom_vec_meso.cpp
action pair_sdpd_taitwater_isothermal.h atom_vec_meso.cpp
action fix_rigid_sph.cpp fix_rigid.cpp
action fix_rigid_sph.h fix_rigid.cpp
action fix_move_sph.cpp atom_vec_meso.cpp
action fix_move_sph.h atom_vec_meso.h
