#!/bin/bash

WGET_VERSION=$(wget --version | grep -oie "wget [0-9][0-9.]*" | head -n 1 | awk '{print $2}')
if [ -z "$WGET_VERSION" ]
then
WGET_VERSION=PARSE_ERROR
fi

WGET_USER_AGENT="wget/$WGET_VERSION/esg/3.0.51-20220308-212225"

safe_download() {

if [ ! -f $(basename $1) ]; then
  download $1
else
    if [ $((`du $(basename $1) | awk '{print $1}'`)) == $2 ]; then
      echo "$(basename $1) already downloaded"
    else
      echo "File exists, but is the incorrect size. Downloading again."
      rm $(basename $1)
      download $1
    fi
fi
}

safe_extract() {
dir_name=`echo $1 | awk '{print substr($1, 1, length($1)-7)}'`
if [ ! -d $dir_name ] ; then
   tar -xvf $1
else
  dir_size=`du -s $dir_name | awk '{print $1}'`
  found=0
  for valid_size in "${@:2}"
  do
    if [[ $valid_size == $dir_size ]]; then
       found=1
    fi
  done
  if [[ $found == 0 ]]; then
      echo "Directory $dir_name exists, but is an unexpected size. Clearing directory and unpacking again."
      rm -rf $dir_name
      tar -xvf $1
  else
      echo "$dir_name already unpacked."
  fi
fi
}

safe_extract_with_dirname() {
dir_name=$1
if [ ! -d $dir_name ] ; then
   tar -xvf $2
else
  dir_size=`du -s $dir_name | awk '{print $1}'`
  found=0
  for valid_size in "${@:3}"
  do
    if [[ $valid_size == $dir_size ]]; then
       found=1
    fi
  done
  if [[ $found == 0 ]]; then
      echo "Directory $dir_name exists, but is an unexpected size. Clearing directory and unpacking again."
      rm -rf $dir_name
      tar -xvf $2
  else
      echo "$dir_name already unpacked."
  fi
fi
}

safe_download_and_extract ()
{
safe_download $1 $2
safe_extract $(basename $1) "${@:3}"
}

download() {
   echo "Downloading '$(basename $1)'..."
   if [[ "$cert" ]]; then
      wget="wget --no-check-certificate -c --user-agent=$WGET_USER_AGENT"
    else
      wget="wget -c --user-agent=$WGET_USER_AGENT"
    fi

    $wget $1 || { failed=1; break; }
}

safe_build() {
if [ ! -f $1 ]; then
  $2
else
    if [ $((`du $1 | awk '{print $1}'`)) == $3 ]; then
      echo "$1 has already been built."
    else
      echo "$1 exists, but is the incorrect size. Compiling again..."
      $2
    fi
fi
}


download_hdf5() {
   echo "Downloading 'HDF5' code..."
   if [[ "$cert" ]]; then
      wget="wget --no-check-certificate -c --user-agent=$WGET_USER_AGENT"
    else
      wget="wget -c --user-agent=$WGET_USER_AGENT"
    fi

    $wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_12_1.tar.gz || { failed=1; break; }
}


download_cdo() {
   echo "Downloading 'Climate Data Operators' code..."
   if [[ "$cert" ]]; then
      wget="wget --no-check-certificate -c --user-agent=$WGET_USER_AGENT"
    else
      wget="wget -c --user-agent=$WGET_USER_AGENT"
    fi

    $wget https://code.mpimet.mpg.de/attachments/download/26823/cdo-2.0.5.tar.gz || { failed=1; break; }
}

build_sqlite()
{
    start_dir=$PWD
    cd sqlite-version-3.38.2
    ./configure --prefix=$1
    make -j 99 && make install
    cd $start_dir
}

build_proj()
{
    start_dir=$PWD
    cd proj-8.2.0
    ./configure --prefix=$1 --without-curl --without-mutex
    make -j 99 && make install
    cd $start_dir
}

build_hdf5() {
    start_dir=$PWD
    cd hdf5-hdf5-1_12_1
    ./configure --prefix=$1
    make -j 99 && make install
    cd $start_dir
}

build_netcdf() {
    start_dir=$PWD
    cd netcdf-c-4.8.1
    ./configure --prefix=$1 --enable-shared --enable-netcdf-4 --with-pic
    make -j 99 && make install
    cd $start_dir
}

build_cdo() {
    start_dir=$PWD
    cd cdo-2.0.5
    ./configure --with-netcdf=$1 --with-hdf5=$1 --with-proj=$1 --prefix=$1
    make -j 99 && make install
    cd $start_dir
}

#
# MAIN
#
start_dir=$PWD

# build cdo
mkdir -p build_cdo
cd build_cdo

mkdir -p install
install_dir=$PWD/install
safe_download_and_extract https://code.mpimet.mpg.de/attachments/download/26823/cdo-2.0.5.tar.gz 11516 647200 647232
safe_download https://github.com/sqlite/sqlite/archive/refs/tags/version-3.38.2.tar.gz 11728
safe_extract_with_dirname sqlite-version-3.38.2 version-3.38.2.tar.gz 107392 107388
safe_download_and_extract https://download.osgeo.org/proj/proj-8.2.0.tar.gz 5744 29112
safe_download https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_12_1.tar.gz 13600
safe_extract_with_dirname hdf5-hdf5-1_12_1 hdf5-1_12_1.tar.gz 129364 221228
safe_download https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.8.1.tar.gz 18516
safe_extract_with_dirname netcdf-c-4.8.1 v4.8.1.tar.gz 110244 132472

safe_build $install_dir/bin/sqlite3 "build_sqlite $install_dir" 6928

export LDFLAGS="-L$install_dir/lib -lsqlite3"
export SQLITE3_LIBS="-L$install_dir/lib -lsqlite3"

export CPPFLAGS=-I/$install_dir/include
export SQLITE3_CFLAGS=-I/$install_dir/include

safe_build $install_dir/bin/proj "build_proj $install_dir" 216
safe_build $install_dir/bin/h5dump "build_hdf5 $install_dir" 264
safe_build $install_dir/bin/ncdump "build_netcdf $install_dir" 116
safe_build $install_dir/bin/cdo "build_cdo $install_dir" 107564
