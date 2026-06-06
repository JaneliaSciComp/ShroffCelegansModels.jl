# Base image is overridden by each namespace's BuildConfig (dockerStrategy.from)
FROM scratch

COPY . /app/ShroffCelegansModels.jl

# Re-precompile only the local package and anything that depends on it.
# External deps are already precompiled in the base image.
RUN julia --threads=auto --project=/app/ShroffCelegansModels.jl/web \
    -e 'using Pkg; Pkg.precompile()'

CMD ["julia", "--version"]
