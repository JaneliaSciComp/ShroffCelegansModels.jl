FROM julia:1.12

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

COPY . /app/ShroffCelegansModels.jl

# Precompile all packages (including local ShroffCelegansModels dep).
# Build takes ~15-30 min; container startup is fast thereafter.
RUN julia --threads=auto --project=/app/ShroffCelegansModels.jl/web \
    -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

CMD ["julia", "--version"]
