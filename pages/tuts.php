<?php
   if ($psub=="") { $psub="summary"; }
   
  function submenu($item,$label)
  {
     global $page; global $psub;
     echo "<a class='submenu ".($item==$psub?"subsel":"submenu")."' href='compose.php?page=".$page."&amp;psub=".$item."'>".$label."</a>";
  }
?>
   <div class="columns">
     <div class="menu" style="margin:10px -10px 0px -10px">
     <?php    
       submenu("summary","summary");
       submenu("analysis","analysis");
       submenu("landmarks","landmarks");
       submenu("sketch-map","sketch-map");  
       submenu("results","results");
       submenu("projecting","projecting");       
       submenu("gnuplot","gnuplot");        
     ?>
     </div>
     <?php         
       include("pages/tuts/".$psub.".html")
     ?>
   </div>
