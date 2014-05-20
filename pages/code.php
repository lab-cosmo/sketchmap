<div class="columns">
   <h2> Sketch-map tools </h2>
   <p> From these pages you can obtain instructions to download and compile
   a set of tools which are useful to analyze high-dimensional data and 
   perform sketch-map non-linear dimensionality reduction. These tools are
   much of a work in progress and far from perfect. You can download
   the <a href='data/sketchmap-0.2.tar.bz2'>source as a tarball</a>, but
   it is recommended that you get the latest development version using 
   <span class="code">git</span>, as explained below. 
   </p>
   <h3> Get the source using git </h3>

   <p> The latest development version of the sketch-map tools is stored in a git repository
   you can access this repository by issuing the following command.
   </p>

   <div class="code">
   git clone git://git.berlios.de/sketchmap
   </div>

   <p>
   This will crate a directory called sketchmap that contains the latest version of the code. 
   For more details on how to keep your distribution up to date visit the <a href="http://git-scm.com/">git website</a>. 
   </p>
   
   <h3> Compilation instructions </h3>

   <p>Compiling sketch-map is straightforward.  Once you have downloaded the tar ball simply issue the following commands:</p>

   <div class="code">
   tar -xvf sketchmap-0.2.tar.bz2
   cd sketchmap-0.2
   make
   </div>

   <p>Once the first of these commands has been issued a new directory called sketchmap-0.1 should be created. This contains 
   the source files for the sketchmap tools and the Makefile that tells your computer how to compile the source. To compile
   the code on your computer it may be necessary to edit the make.in file. For example if you don't want to use the 
   g++ compiler then you should change the compiler inside the make.in file. To be clear though oftentimes you can compile the
   without making any changes to make.in.  In our experience the only thing that can sometimes need changing is the location
   of the lapack library.  If you are compiling the code on a unix machine you should set this using:</p>
   </p>
   <div class="code">
   LLAPACK=-llapack
   </div>    

   <p>
   On a mac though this should be set equal to:
   </p>
   <div class="code">
   LLAPACK=-framework veclib
   </div>

   <h3> Installing GISMO </h3>

   <p>When you download the sketch-map repository you will find a directory inside the tar ball called gismo.  This directory contains the
   code that controls the graphical user interface.  To add this code to an existing <a href="http://www.ks.uiuc.edu/Research/vmd/">vmd</a> installation:
   </p>
   
   <ul>
   <li>
   Find where the vmd libraries are installed on your computer.  By default these will be in /usr/local/lib/vmd.
   In the vmd directory there should be a subdirectory called scripts.  If there is not already one there you will
   need to make a new directory inside vmd/scripts called init.d. <br/> <br/>
   </li>
   <li>
   Copy the contents of the gismo directory (all the *.tcl files) from the unpacked sketch-map tar ball to your new init.d directory inside vmd.
   GISMO will now load automatically when you load vmd <br/> <br/>
   </li>
   <li>
   To take full advantage of GISMO you should also download <a href="http://www.plumed-code.org">plumed</a> and compile the sum_hills and driver utilities.
   To run these programs directly from GISMO you must set an environmental variable, plumedir, inside your .bashrc file so that GISMO can find the top directory
   of your plumed distribution. <br/> <br/>
   </li>
   <li>
   To use the sketch-map tools in GISMO you must set a further environmental variable, smapdir, inside the .bashrc file.
   This directory tells plumed where you have unpacked the sketch-map tar ball and thus where the sketch-map executables are located.
   </li>
   </ul>
</div>
