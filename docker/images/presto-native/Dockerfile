FROM ubuntu:24.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ca-certificates git gcc g++ gfortran make meson ninja-build pkg-config libfftw3-dev libcfitsio-dev libglib2.0-dev libpng-dev libx11-dev pgplot5 time python3 python3-numpy python3-scipy && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/scottransom/presto.git /presto && git -C /presto checkout 7f66332b03c52a55603308f755958f030da0c964
ENV PRESTO=/presto
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends libgsl-dev liberfa-dev && rm -rf /var/lib/apt/lists/*
RUN cd /presto && meson setup build --prefix=/opt/presto -Dc_args=-mcpu=native && meson compile -C build && meson install -C build
ENV LD_LIBRARY_PATH=/opt/presto/lib:/opt/presto/lib/aarch64-linux-gnu PATH=/opt/presto/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin LIBRARY_PATH=/opt/presto/lib:/opt/presto/lib/aarch64-linux-gnu
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-pip python3-dev python3-matplotlib && rm -rf /var/lib/apt/lists/*
RUN cd /presto/python && PKG_CONFIG_PATH=/opt/presto/lib/aarch64-linux-gnu/pkgconfig pip install --break-system-packages --ignore-installed --no-cache-dir .
RUN ln -sf /usr/bin/python3 /usr/local/bin/python
