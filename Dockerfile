FROM python:3.11-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:${PATH}"
# Set library paths for Debian architecture detection
ENV LD_LIBRARY_PATH="/usr/local/lib:/usr/lib/$(dpkg-architecture -qDEB_HOST_MULTIARCH)"
ENV PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:/usr/lib/$(dpkg-architecture -qDEB_HOST_MULTIARCH)/pkgconfig"

# Install build dependencies for Debian slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    make \
    m4 \
    wget \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    cmake \
    git \
    csh \
    tcsh \
    gfortran \
    libhdf5-dev \
    libnetcdf-dev \
    libnetcdf-c++4-dev \
    netcdf-bin \
    pkg-config \
    dpkg-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

# Build and install netcdf-fortran from source
WORKDIR /tmp
RUN wget https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz && \
    tar -xzf v4.6.1.tar.gz && \
    cd netcdf-fortran-4.6.1 && \
    ./configure --prefix=/usr/local --enable-shared --disable-static && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    cd / && rm -rf /tmp/netcdf-fortran-4.6.1*

# Set HYCOM environment variables for Debian multi-arch
ENV NETCDF_DIR=/usr/local \
    NETCDF_INC=/usr/local/include \
    NETCDF_LIB=/usr/local/lib \
    HDF5_DIR=/usr \
    HDF5_INC=/usr/include

# Set HDF5_LIB dynamically based on architecture
RUN export HDF5_LIB="/usr/lib/$(dpkg-architecture -qDEB_HOST_MULTIARCH)" && \
    echo "export HDF5_LIB=${HDF5_LIB}" >> /etc/environment

# Verify NetCDF installation
RUN nc-config --all && nf-config --all

# Copy HYCOM-tools source
WORKDIR /usr/local
# COPY . HYCOM-tools/
RUN git clone https://github.com/HYCOM/HYCOM-tools.git /usr/local/HYCOM-tools
WORKDIR /usr/local/HYCOM-tools

# Fix the config file for gfortran compilation with architecture-aware flags
RUN cp config/intelGF_setup config/intelGF_setup.bak
RUN EXTRANCDF_FLAGS="$(nf-config --flibs) $(nc-config --libs)" && \
    # Set architecture-specific flags
    ARCH_FLAGS="" && \
    if [ "$(uname -m)" = "x86_64" ]; then \
      ARCH_FLAGS="-m64"; \
    fi && \
    echo "#" > config/intelGF_setup && \
    echo "# Linux gfortran configuration for HYCOM (multi-arch)" >> config/intelGF_setup && \
    echo "#" >> config/intelGF_setup && \
    echo "FC            = gfortran" >> config/intelGF_setup && \
    echo "CC            = gcc" >> config/intelGF_setup && \
    echo "LD            = \$(FC)" >> config/intelGF_setup && \
    echo "ARCH          = intelGF" >> config/intelGF_setup && \
    echo "TYPE          = gfortran" >> config/intelGF_setup && \
    echo "" >> config/intelGF_setup && \
    echo "CPPFLAGS      = -DIA32 -DREAL4 -DENDIAN_IO -DTIMER -DRELO" >> config/intelGF_setup && \
    echo "FCFFLAGS      = -fPIC ${ARCH_FLAGS} -fno-second-underscore -O -std=legacy -fallow-argument-mismatch" >> config/intelGF_setup && \
    echo "FFLAGS        = \$(FCFFLAGS) \$(CPPFLAGS) -I/usr/local/include" >> config/intelGF_setup && \
    echo "LDFLAGS       = \$(FFLAGS)" >> config/intelGF_setup && \
    echo "CFLAGS        = -O -I/usr/local/include" >> config/intelGF_setup && \
    echo "" >> config/intelGF_setup && \
    echo "EXTRANCDF     = ${EXTRANCDF_FLAGS}" >> config/intelGF_setup

# Fix Make_all.src to use gfortran instead of amdIF
RUN sed -i 's/setenv ARCH amdIF/setenv ARCH intelGF/' Make_all.src

# Create proper Make_ncdf.src configuration  
RUN EXTRANCDF_FLAGS="$(nf-config --flibs) $(nc-config --libs)" && \
    echo "# NetCDF configuration for gfortran" > Make_ncdf.src && \
    echo "setenv ARCH intelGF" >> Make_ncdf.src && \
    echo "setenv NCDFC /usr" >> Make_ncdf.src && \
    echo "setenv NCDF /usr/local" >> Make_ncdf.src && \
    echo "setenv EXTRANCDF \"${EXTRANCDF_FLAGS}\"" >> Make_ncdf.src

# Ensure Make_all.csh uses LinuxGF (GNU Fortran) rather than LinuxAIF inside container
RUN sed -E -i 's/^#?[ \t]*setenv OS LinuxGF$/setenv OS LinuxGF/' bin/Make_all.csh && \
    sed -E -i 's/^([ \t]*)setenv OS LinuxAIF$/# \1setenv OS LinuxAIF/' bin/Make_all.csh

    # Ensure Make_ncdf.csh uses LinuxGF (GNU Fortran) rather than LinuxAIF inside container
RUN sed -E -i 's/^#?[ \t]*setenv OS LinuxGF$/setenv OS LinuxGF/' bin/Make_ncdf.csh && \
    sed -E -i 's/^([ \t]*)setenv OS LinuxAIF$/# \1setenv OS LinuxAIF/' bin/Make_ncdf.csh

# Build base tools first
RUN csh Make_all.csh 2>&1 | tee make_all.log

# Build NetCDF-enabled tools with explicit configuration
RUN csh Make_ncdf.csh 2>&1 | tee make_ncdf.log || true

# Install Python dependencies (pip is already available in python:3.11-slim)
RUN pip install --no-cache-dir xarray netcdf4 dask numpy awscli

# Set working directory
WORKDIR /data

# Create symbolic links for easier access
RUN ln -sf /usr/local/HYCOM-tools/archive/src/archv2ncdf2d /usr/local/bin/ && \
    ln -sf /usr/local/HYCOM-tools/archive/src/archv2ncdf3z /usr/local/bin/

CMD ["tail", "-f", "/dev/null"]