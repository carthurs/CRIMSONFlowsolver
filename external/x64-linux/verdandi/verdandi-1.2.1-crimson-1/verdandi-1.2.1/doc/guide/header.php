<?php echo '<?xml version="1.0"  encoding="iso-8859-1"?'.'>' ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<?php $root='..';?>

<head>
<meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
<meta name="keywords" content="Verdandi, library, C++, Python, data,
assimilation, data assimilation, estimation, parameter estimation, filter,
Kalman, sequential, variational, uncertainty, aggregation, implementation,
open source"/>
<meta name="description" content="Verdandi is a generic C++ library for data
assimilation. It aims at providing methods and tools for data assimilation. It
is designed to be relevant to a large class of problems involving
high-dimensional numerical models."/>
<title>Verdandi user's guide</title>
<link rel="stylesheet" type="text/css" href="<?php echo $root?>/content.css"/>
<link rel="stylesheet" href="tabs.css" type="text/css"/>
<link rel="stylesheet" href="guide.css" type="text/css"/>
<?php if (file_exists($root.'/prettify.js'))
  echo '<script type="text/javascript" src="'.$root.'/prettify.js"></script>';
else if (file_exists('prettify.js'))
  echo '<script type="text/javascript" src="prettify.js"></script>'; ?>
</head>

<body onload="prettyPrint()">

<div class="page">

<?php if (file_exists($root.'/header.php'))
      include $root.'/header.php'; ?>

<div class="doc">

<?php function HL($file_, $section_, $string_)
{
if ($file_ == $section_)
  echo '<em>'.$string_.' </em>';
else
  echo '<a href="'.$section_.'.php">'.$string_.'</a>';
}; ?>

<?php $file=basename($_SERVER['REQUEST_URI'], ".php"); $file = explode(".", $file); $file = $file[0];?>

<div class="nav">

<ul>
<li class="jelly"> <b>USER'S GUIDE</b> </li>
<li class="jelly"> <?php HL($file, "index", "Introduction");?>  </li>

<li class="jelly"> <?php HL($file, "getting_started", "Getting Started");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "getting_started"
or basename($_SERVER['REQUEST_URI'], ".php") == "installation"
or basename($_SERVER['REQUEST_URI'], ".php") == "example_programs"
or basename($_SERVER['REQUEST_URI'], ".php") == "notation"
or basename($_SERVER['REQUEST_URI'], ".php") == "overview")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "installation", "Installation");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "example_programs", "Example Programs");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "notation", "Notation");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "overview", "Overview");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "assimilation_methods", "Assimilation Methods");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "assimilation_methods"
or basename($_SERVER['REQUEST_URI'], ".php") == "optimal_interpolation"
or basename($_SERVER['REQUEST_URI'], ".php") == "extended_kalman_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "reduced_order_extended_kalman_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "unscented_kalman_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "reduced_order_unscented_kalman_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "reduced_minimax_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "monte_carlo"
or basename($_SERVER['REQUEST_URI'], ".php") == "ensemble_kalman_filter"
or basename($_SERVER['REQUEST_URI'], ".php") == "four_dimensional_variational")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "optimal_interpolation", "Optimal Interpolation");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "extended_kalman_filter", "Extended Kalman Filter");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "reduced_order_extended_kalman_filter", "Reduced Order Extended Kalman Filter");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "unscented_kalman_filter", "Unscented Kalman Filter");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "reduced_order_unscented_kalman_filter", "Reduced Order Unscented Kalman Filter");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "reduced_minimax_filter", "Reduced Minimax Filter");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "monte_carlo", "Monte Carlo");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "ensemble_kalman_filter", "Ensemble Kalman Filter");
  echo '</li>';
  HL($file, "four_dimensional_variational", "Four Dimensional Variational");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "models", "Models");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "models"
or basename($_SERVER['REQUEST_URI'], ".php") == "quadratic_model"
or basename($_SERVER['REQUEST_URI'], ".php") == "shallow_water_model"
or basename($_SERVER['REQUEST_URI'], ".php") == "clamped_bar_model"
or basename($_SERVER['REQUEST_URI'], ".php") == "lorenz_model")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "quadratic_model", "Quadratic Model");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "shallow_water_model", "Shallow-water");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "clamped_bar_model", "Clamped Bar");
  echo '</li>';
  HL($file, "lorenz_model", "Lorenz model");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "observations", "Observations");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "observations"
or basename($_SERVER['REQUEST_URI'], ".php") == "linear_observation_manager"
or basename($_SERVER['REQUEST_URI'], ".php") == "grid_to_network_observation_manager"
or basename($_SERVER['REQUEST_URI'], ".php") == "observation_aggregator")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "linear_observation_manager", "Linear Observation Manager");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "grid_to_network_observation_manager", "Grid To Network Observation Manager");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "observation_aggregator", "Observation Aggregator");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "tools", "Tools");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "tools"
or basename($_SERVER['REQUEST_URI'], ".php") == "perturbation_manager")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "perturbation_manager", "Perturbation Manager"); 
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "plugging_in_verdandi", "Plugging in Verdandi");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "plugging_in_verdandi"
or basename($_SERVER['REQUEST_URI'], ".php") == "plugging_model"
or basename($_SERVER['REQUEST_URI'], ".php") == "plugging_observation")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "plugging_model", "Model");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "plugging_observation", "Observations");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "using_verdandi", "Using Verdandi");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "using_verdandi"
or basename($_SERVER['REQUEST_URI'], ".php") == "tips"
or basename($_SERVER['REQUEST_URI'], ".php") == "configuration_files"
or basename($_SERVER['REQUEST_URI'], ".php") == "scons")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "configuration_files", "Lua Configuration Files");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "scons", "SCons");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "tips", "Tips");
  echo '</li> </ul>';
} ?>
</li>

<li class="jelly"> <?php HL($file, "debugging", "Debugging");?>  </li>

<li class="jelly"> <?php HL($file, "python", "Python");?> </li>
<li class="jelly"> <b>API REFERENCE</b> </li>
<li class="jelly"> <?php HL($file, "annotated", "Classes");?>
<ul class="navsubul"> <li class="jelly"> <?php HL($file, "annotated", "Class List");?> </li> 
<li class="jelly"> <?php HL($file, "hierarchy", "Class Hierarchy");?> </li>
<li class="jelly"> <?php HL($file, "functions", "Class Members");?>
</li> </ul> </li>
<li class="jelly"> <?php HL($file, "namespacemembers", "Functions");?> </li>
<li class="jelly"> Search for <form action="search.php" method="get">
    <input class="search" type="text" name="query" value="" size="20" accesskey="s"/>
  </form>
</li>
<!-- <li class="jelly"> <?php HL($file, "faq", "F.A.Q.");?> </li>-->
<li class="jelly"> <a
href="mailto:verdandi-help@lists.gforge.inria.fr"
style="color:black">Support</a></li>
</ul>

</div>

<div class="doxygen">
<div>
