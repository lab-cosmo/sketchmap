<div class="columns">
<!--
  <div class="lcolumn">
-->
    <img src="images/gismo_logo.png"  width=400 alt="A galary of gismo images" class="floatright"/>
    <p>
     GISMO is a tcl package that  can help one to perform dimensionality reduction calculations and to visualize their results.  
     It is designed to incorporate into <a href="http://www.ks.uiuc.edu/Research/vmd/">vmd</a> - in fact it <b>will not</b> work independently of vmd.  vmd is a package for visualizing atomic trajectories.  
     Each frame in a trajectory is a high dimensionality vector that describes the position of many atoms.  Using GISMO one can create a map that 
     connects every frame to to a point in a low-dimensionality space.  These low-dimensional vectors thus constitute a simpler representation of 
     the trajectory.  GISMO then allows one to visualize this mapping through a simple point and click interface.
    </p>
    <p>
     In atomistic simulations dimensionality reduction is routinely performed by dreaming up chemically motivated collective 
     variables (CVs).  These could be distances between particularly interesting atoms, angles, torsions or some more complicated 
     function.  GISMO can calculate many CVs as it contains an interface to the <a href="http://www.plumed-code.org"> plumed </a> package.
    </p>
    <p>
     Manifold learning algorithms provide an alternative to CVs and are much less reliant on chemically intuition.  GISMO contains 
     a graphical user interface for doing this sort of analysis with sketch-map and for visualizing the results.
    </p>
     GISMO has been written by Gareth Tribello and is still far from complete.  Users are kindly requested to report any bugs they find to
     gareth dot tribello at phys dot chem dot ethz dot ch.
    </p>
<!--
  </div>
  <div class="rcolumn">
     <br> <br>
-->
<!--
  </div>
</div>
<div class="columns">
-->

<h3> Getting started </h3>
<!--
</div>
<div class="columns">
  <div class="lcolumn">
-->
    <img src="images/gismo_window.png"  width=300 alt="The main gismo window" class="floatleft"/>
<!--
  </div>
  <div class="rcolumn"> 
    <br>  <br>
-->
    <p>
    You can open GISMO by clicking on Extensions/GISMO in the vmd main window.  The window shown in the left column will then appear.  The controls in the upper
    menu bar provide allow you to load and save data (File), to calculate properties of your trajectory (Calculate) or to adjust the appearance of your plot (Appearance). 
    Meanwhile the bottom menu bar controls what data is plotted in the main window.  GISMO is able to plot two kinds of data:
    </p>

    <ul>
    <li> Collective variable (CVs) data </li>
    <li> One or two dimensional free energy surfaces </li>
    </ul>

    <p>
    These quantities can be shown simultaneously on the same plot or displayed separately by manipulating the controls in the bottom pannel.
    </p>
<!--
  </div>
</div>
<div class="columns">
-->
<h4> Collective variable mode </h4>

<p>
Once you have selected quantities to display on the x-axis and y-axis every point in the trajectory for the top molecule is represented with 
a single square in the 2d plane.  The sizes and colors of these squares can then be used to display the values of further collective variables.  
One square is colored black.  This colored square is the one corresponding to trajectory frame from the top molecule that is displayed in the 
main vmd window.  As you move through the trajectory this colored square will be updated in accordance with the trajectory.  Furthermore, if you
click on any one of the squares displayed in the two dimensional plane the corresponding trajectory frame will be shown in the main vmd window.
</p>

<p>
If there is some particularly intereting portion of CV space that you want to look at then drag a box around that area whilst holding down the 
shift key and the left mouse button.  GISMO will then zoom in on the selected area.
</p>

<p>
Alternatively, if you wish to select the frames in a particular portion of CV space, then you can do this by draging a box around the area of interest
whilst holding down the left mouse button and the control key.  When you release the mouse button all the frames in that area will be colored red.  
You can save these frames in a new molecule by selecting Calculate/Save selection.    
</p>

<h4> Free energy mode </h4>

<p>
You can load a free energy surface or alternatively calculate a free energy surface using plumed's sum_hills utility.  The gui keeps a record of all
the free energy surfaces you load/calculate.  The free energy surface shown is controlled by adjusting the CVs plotted and through the fes data menu
in the bottom pannel.  Furthermore, GISMO provides simple tools for visualizing how the free energy surface change as a function of simulation time 
(the fes frame counter).  This is very useful when it comes to examining the convergence of free energies.
</p>

<p> 
You can integrate the free energy in a basin on your free energy surface (or a free energy difference between two basins) using a simple procedure.
Push and hold the left mouse button while holding down the control key to create a square containing the basin whose free energy you want to integrate.
</p>

<p>
Please be aware that the free energy functions are the part of GISMO where there are most likely to be bugs.  Also the routines can be rather slow especially
when you are plotting two dimensional free energy surfaces.
</p> 

</div>
