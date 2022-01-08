# lgnrmfit
A utility that can fit the lognormal distribution to a number of data

# How to build
lgnrlm depends on the [dlib](https://github.com/davisking/dlib).</br>
Before building __lgnrmfit__ clone dlib.
Let's suppose that we want to clone dlib to the `dlib_dir` directory
```
cd dlib_dir
git clone https://github.com/davisking/dlib.git
```
Then the `DLIB_PATH` is `dlib_dir/dlib/dlib` which is the folder with the source code of dlib.


The following steps download and build the repo
```
git clone https://github.com/giorgk/lgnrmfit.git
cd lgnrmfit
mkdir build
cd build
cmake -DDLIB_PATH=/path/to/dlib/dlib ..
make
```

# How to use
:construction:
