<!-- Required extensions: pymdownx.betterem, pymdownx.tilde, pymdownx.emoji, pymdownx.tasklist, pymdownx.superfences -->
<!-- The first line is for ReText editor -->

# Minimal guide to compiling LKH TSP solver for VeRyPy

The guide assumes Unbutu/linux. For Windows use [MSYS2](https://www.msys2.org/) or run VeRyPy entirely inside [WSL](https://docs.microsoft.com/en-us/windows/wsl) environment.

First, create the path and folders for the TSP solvers.
One can use whatever path, but using the predefined path removes the need to modify configuration files. Then, download the LKH source code using `wget` and uncompress the archive.

```console
$ mkdir -p ~/Projects/Research/TSP
$ cd ~/Projects/Research/TSP 
$ wget http://akira.ruc.dk/~keld/research/LKH/LKH-2.0.9.tgz
$ tar -xzvf LKH-2.0.9.tgz 
$ cd LKH-2.0.9/
```

For the next step you need to install C compiler (cc) and make.
On Ubuntu and WSL these tools are provided by the build-essential metapackage.

```console
$ sudo apt install build-essential
```

For Windows / MSYS2 there is a similar package called base-devel.
```console
$ pacman -S base-devel
```

When you have a C development environment installed, you can proceed and compile the LKH.

```console
$ make
```

Test the still hot LKH executable

```console
$ ./LKH pr2392.par pr2392.tsp
```

If everything was succesfull, you get some processing, output, and eventually the result.

```console
Successes/Runs = 10/10
Cost.min = 378032, Cost.avg = 378032.00, Cost.max = 378032
Gap.min = 0.0000%, Gap.avg = 0.0000%, Gap.max = 0.0000%
Trials.min = 1, Trials.avg = 5.8, Trials.max = 29
Time.min = 0.21 sec., Time.avg = 0.47 sec., Time.max = 0.98 sec.
```

Note that similar process also works for compiling ACOTSP. It is used as a very fast TSP solver for CMT-2P.
The ACOTSP wrapper in VeRyPy expects the custom version from [https://github.com/yorak/ACOTSP](https://github.com/yorak/ACOTSP),
which allows disabling the ant systems and only the local search to be used.
 (which can be 

If you used different path or version of the LKH, please be sure to modify the TSP executable paths in `config.py` of VeRyPy.
