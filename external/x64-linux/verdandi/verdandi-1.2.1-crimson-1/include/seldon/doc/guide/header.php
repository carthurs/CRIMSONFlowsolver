<?php echo '<?xml version="1.0"  encoding="iso-8859-1"?'.'>' ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<?php $root='..';?>

<head>
<meta http-equiv="Content-Type" content="text/html;charset=UTF-8">
<title>Seldon user's guide</title>
<link rel="stylesheet" type="text/css" href="<?php echo $root?>/content.css">
<link rel="stylesheet" href="tabs.css" type="text/css">
<link rel="stylesheet" href="guide.css" type="text/css">
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
<li class="jelly"> <?php HL($file, "installation", "Installation");?> </li>
<li class="jelly"> <?php HL($file, "overview", "Overview");?> </li>
<li class="jelly"> <?php HL($file, "vectors", "Vectors");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "vectors"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_vector"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_sparse_vector"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_vector")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "class_vector", "Dense Vectors");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "class_sparse_vector", "Sparse Vectors");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_vector", "Functions");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "matrices", "Matrices");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "matrices"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_matrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "class_sparse_matrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_matrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "submatrix"
or basename($_SERVER['REQUEST_URI'], ".php") == "matrix_miscellaneous")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "class_matrix", "Dense Matrices");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "class_sparse_matrix", "Sparse Matrices");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_matrix", "Functions");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "submatrix", "Sub-Matrices");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "matrix_miscellaneous", "Miscellaneous");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "other_structures", "Other Structures");?>
<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "other_structures"
or basename($_SERVER['REQUEST_URI'], ".php") == "vector2"
or basename($_SERVER['REQUEST_URI'], ".php") == "array3d")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "vector2", "Vector2");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "array3d", "3D&nbsp;Array");
  echo '</li> </ul>';
} ?>
</li>

<li class="jelly"> <?php HL($file, "allocators", "Allocators");?>  </li>
<li class="jelly"> <?php HL($file, "exceptions", "Exceptions");?>  </li>
<li class="jelly"> <?php HL($file, "computations", "Computations");?>

<?php if (basename($_SERVER['REQUEST_URI'], ".php") == "computations"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_blas"
or basename($_SERVER['REQUEST_URI'], ".php") == "functions_lapack"
or basename($_SERVER['REQUEST_URI'], ".php") == "direct"
or basename($_SERVER['REQUEST_URI'], ".php") == "eigenvalue"
or basename($_SERVER['REQUEST_URI'], ".php") == "iterative")
{
  echo '<ul class="navsubul"> <li class="jelly">';
  HL($file, "functions_blas", "Blas");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "functions_lapack", "Lapack");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "direct", "Direct Solvers");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "iterative", "Iterative Solvers");
  echo '</li>';
  echo '<li class="jelly">';
  HL($file, "eigenvalue", "Eigenvalue Solvers");
  echo '</li> </ul>';
} ?>

</li>

<li class="jelly"> <?php HL($file, "python", "Python Interface");?> </li>
<li class="jelly"> <?php HL($file, "glossary", "Index");?> </li>
<li class="jelly"> <b>API REFERENCE</b> </li>
<li class="jelly"> <?php HL($file, "annotated", "Classes");?>
<ul class="navsubul"> <li class="jelly"> <?php HL($file, "annotated", "Class List");?> </li> 
<li class="jelly"> <?php HL($file, "hierarchy", "Class Hierarchy");?> </li>
<li class="jelly"> <?php HL($file, "functions", "Class Members");?>
</li> </ul> </li>
<li class="jelly"> <?php HL($file, "namespacemembers", "Functions");?> </li>
<li class="jelly"> Search for <form action="search.php" method="get">
    <input class="search" type="text" name="query" value="" size="20" accesskey="s">
  </form>
</li>
<!-- <li class="jelly"> <?php HL($file, "faq", "F.A.Q.");?> </li>-->
<li class="jelly"> <a
href="mailto:seldon-help@lists.sourceforge.net"
style="color:black">Support</a></li>
</ul>

</div>

<div class="doxygen">
