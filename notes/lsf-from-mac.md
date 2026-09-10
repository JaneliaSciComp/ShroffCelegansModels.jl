# Using LSF from a Mac

The Janelia LSF cluster (`bsub`, `bjobs`, `bhist`) is only available on Janelia
Linux login nodes. A Mac cannot run these commands directly; use SSH to proxy
them through a login node.

## One-off commands via SSH

```bash
# Check a job
ssh login1.int.janelia.org bjobs 12345

# Submit a job script that is already on the shared filesystem
ssh login1.int.janelia.org "cd /groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl && bsub < scripts/submit_movie_glmakie.bsub"

# Job history
ssh login1.int.janelia.org bhist -l 12345
```

## Interactive login shell

```bash
ssh login1.int.janelia.org
# then run bsub / bjobs / bhist as normal
```

## Submitting a job that reads files on the shared filesystem

The cluster compute nodes and login nodes share the NFS filesystem at
`/groups/…`. If your input files (e.g. `*.h5`) are already there you can
submit directly from the Mac via SSH (see one-off command pattern above).

## Accessing `/groups/…` from Mac

The `/groups/` NFS filesystem is accessible on Mac via **SMB mount**:

- In Finder: `Go → Connect to Server…` → `smb://groups.hhmi.org/scicompsoft`
- Or from Terminal: `open smb://groups.hhmi.org/scicompsoft`

Once mounted it appears as `/Volumes/scicompsoft/`. This lets you read log
files and copy input/output files without SSH file transfer commands.

## Reading log files from Mac

After mounting via SMB, log files written to `/groups/scicompsoft/…` are
directly readable at the corresponding `/Volumes/scicompsoft/…` path.

Alternatively, tail them over SSH:

```bash
ssh login1.int.janelia.org \
  tail -f /groups/scicompsoft/home/kittisopikulm/src/ShroffCelegansModels.jl/movies/job_12345.log
```

## Network access

LSF, the login nodes, and the SMB mounts are only reachable inside the Janelia
campus network or over the **Janelia VPN**. Ensure the VPN is connected before
any of the above.
