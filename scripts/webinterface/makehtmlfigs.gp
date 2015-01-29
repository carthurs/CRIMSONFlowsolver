is_missing(x)=int(system("ismissing.sh ".x))

set terminal svg size 1024,768 enhanced font "Helvetica,20"
set output 'PLV.svg'
if (! is_missing("PLV.dat") ) plot './PLV.dat' u ($2/1333.3) w l

set output 'QMV.svg'
if (! is_missing("QMV.dat") ) plot './QMV.dat' u 2 w l

set output 'VLV.svg'
if (! is_missing("VLV.dat") ) plot './VLV.dat' u 2 w l

set output 'ELV.svg'
if (! is_missing("ELV.dat") ) plot './ELV.dat' u 2 w l

set output 'Pstab.svg'
if (! is_missing("Pstab.dat") ) plot './Pstab.dat' u ($2/1333.3) w l

set output 'Pstab.svg'
if (! is_missing("PStab.dat") ) plot './PStab.dat' u ($2/1333.3) w l

set output 'Paorta.svg'
if (! is_missing("Paorta.dat") ) plot './Paorta.dat' u ($2/1333.3) w l

set output 'Qaorta.svg'
if (! is_missing("Qaorta.dat") ) plot './Qaorta.dat' u 2 w l

set output 'sysAhist.svg'
if (! is_missing("sysAhist.dat") ) plot './sysAhist.dat' u 2 w l

set output 'MVO2_history.svg'
if (! is_missing("MVO2_history.dat") ) plot './MVO2_history.dat' u 1 w l

set output 'r_d_history1.svg'
if (! is_missing("r_d_history.dat") ) plot './r_d_history.dat' u 1 w l

set output 'r_p_history1.svg'
if (! is_missing("r_p_history.dat") ) plot './r_p_history.dat' u 1 w l

set output 'meanO2Discrepancy_history1.svg'
if (! is_missing("meanO2Discrepancy_history.dat") ) plot './meanO2Discrepancy_history.dat' u 1 w l


if (! is_missing("PressHist.dat"))
{
 do for [i=2:PressFlowHistSurfaces] {
   title = sprintf("outletPressure%d.svg",i-1)
   set output title
   plot './PressHist.gplot' u i w l

   title = sprintf("outletFlow%d.svg",i-1)
   set output title
   plot './FlowHist.gplot' u i w l   
  }
}
