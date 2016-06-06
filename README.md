This repo is for voice activity detection algorithm (VAD). This code is based on this paper:

Ramırez, Javier, José C. Segura, Carmen Benıtez, Angel De La Torre, and Antonio Rubio. 
   "Efficient voice activity detection algorithms using long-term speech information." 
                       
And the Matlab code version (Thanks for the isrish making such a beautiful program):

https://github.com/isrish/VAD-LTSD


This cpp file take input audio data from txt file and write results to 'example*.txt' file. At the end, I use gnuplot to plot the result and it works exactly the same as the Matlab Code from isrish.

# Compile and run
You need <a href="http://www.fftw.org"> FFTW </a> and <a href="http://http://eigen.tuxfamily.org/"> Eigen3 </a> libraries.
```
$ mkdir build & cd build
$ cmake .. & make
$ ../bin/./vad
```
## Visualizing the result
```
$ gnuplot
$ > set multiplot
$ > plot "sound.txt" with 1:2 lines
$ > filename(n) = sprintf("example_%d.txt", n)
$ > plot for [i=1:10] filename(i) using 1:2 with lines
```
