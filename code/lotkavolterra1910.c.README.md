# Installation

To compile loktavolterra1910.c, one can install the sundials
dependency using spack:

```console
$ cd
$ git clone https://github.com/spack/spack.git
$ source ~/spack/share/spack/setup-env.sh
$ cd -
$ spack env activate .
$ spack concretize -f
$ spack install
```

Now you can compile the loktavolterra1910.c file using sundials with:

```console
$ spack load sundials
$ mkdir -p ../build-lv
$ cd ../build-lv
$ cmake ../code
$ make
```


# Troubleshooting

The installation may hang from using python supplied by homebrew:

```console
$ spack install
==> Installing environment /Users/pnanda/immunology/papers/calipro-vs-abc/code
==> Installing pkgconf-1.8.0-bunqy7qohf3b6dqlzmhzjewnbhv4exlb
==> No binary for pkgconf-1.8.0-bunqy7qohf3b6dqlzmhzjewnbhv4exlb found: installing from source
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/usr/local/Cellar/python@3.10/3.10.8/Frameworks/Python.framework/Versions/3.10/lib/python3.10/multiprocessing/spawn.py", line 116, in spawn_main
    exitcode = _main(fd, parent_sentinel)
  File "/usr/local/Cellar/python@3.10/3.10.8/Frameworks/Python.framework/Versions/3.10/lib/python3.10/multiprocessing/spawn.py", line 126, in _main
    self = reduction.pickle.load(from_parent)
AttributeError: Can't get attribute 'Mark' on <module 'ruamel.yaml.error' from '/Users/pnanda/Library/Python/3.10/lib/python/site-packages/ruamel/yaml/error.py'>
^C==> Error: Failed to install pkgconf due to KeyboardInterrupt: 

==> Error: Keyboard interrupt.
```

In that case, set an environmental to tell spack to use the system
python instead of a non-system version installed from homebrew, etc:

```console
$ which python3
/usr/local/bin/python3
$ SPACK_PYTHON=/usr/bin/python3
```

Then rerun the concretize step and retry the install.
