# HYCOM-tools

This is the HYCOM-tools repository.

# Where to find the information

All information about the HYCOM-tools (installation, compilation) are described in the [HYCOM-tools Wiki](https://github.com/HYCOM/HYCOM-tools/wiki) page.

# Using HYCOM-Tools Docker Container

The [`Dockerfile`](./Dockerfile) provides a platform-independent (supporting ARM and x86 architectures) container runtime with compiled HYCOM-Tools. The resulting image is based on the Debian linux distro. It can be built using the `docker compose` CLI and the included [`docker-compose`](./docker-compose.yaml) file.

## Docker Quickstart
Run these commands from within the repository directory.

Build the docker image
```bash
docker compose build
```

Run the docker image as a container. The `-d` flag runs the container in the background.
```bash
docker compose up -d
```

Stopping the image: `docker compose down`.

Start an interactive bash session on the image
```bash
docker compose exec -ti hycom bash
```

Example Usage: Run `archv2ncdf3z` conversion tool to process files in the local [data](./data) directory (which is mounted to `/data` in the container).
```bash
docker compose exec -ti hycom bash -c "archv2ncdf3z < /data/archv2ncdf3z.IN"
```

## Full Example
Disk Warning: Running the example script will trigger file downloads totaling ~6GB, requiring, once unzipped, ~18GB of disk space.

Run the [rtofs conversion example script](docs/rtofs_conversion.sh) from within the container to test the tooling for converting full-volume RTOFS binaries to netcdf. You can copy the script to the `./data` folder to access it on the running container.
