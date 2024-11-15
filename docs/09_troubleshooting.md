<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).

When using singularity the following error may occur:

```
RuntimeError: cannot cache function 'rdist': no locator available for file '/home/epi2melabs/...'
```

If you receive this error, you may need to set the numba cache directory to a location that is writable by the singularity 
container. To do this, add the following contents to a file (`numba.config` for example) and use `-c numba.config` in the workflow command
to apply it.

```
profiles {
    singularity {
        singularity {
            enabled = true
            autoMounts = true
            runOptions = '--writable-tmpfs'
        }
    }
}
env {
    NUMBA_CACHE_DIR = "${launchDir}/numba_cache"
}
```
