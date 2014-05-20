<?php
   if ($psub=="") { $psub="tuts0"; }
   
  function submenu($item,$label)
  {
     global $page; global $psub;
     echo "<a class='submenu ".($item==$psub?"subsel":"submenu")."' href='compose.php?page=".$page."&amp;psub=".$item."'>".$label."</a>";
  }
?>
   <div class="columns">
     <div class="menu" style="margin:10px -10px 0px -10px">
     <?php    
       submenu("intro","overview");
       submenu("plumed","plumed"); 
       submenu("dimlandmark","dimlandmark");
       submenu("dimred","dimred");
       submenu("dimproj","dimproj");      
     ?>
     </div>
     <?php         
       include("pages/manu/".$psub.".html")
     ?>
   </div>
