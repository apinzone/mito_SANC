#include "subcell.hpp"
#include "mito.hpp" 

void CSubcell::compute_J_ATP_D(void)
{
  // Corners (no-flux boundary: reflected neighbor used in place of the
  // off-grid neighbor, i.e. the same interior neighbor counted twice)
  J_ATP_D[0+0*nx+0*(nx*ny)]=(ATP_cyto[0+0*nx+0*(nx*ny)+1]+ATP_cyto[0+0*nx+0*(nx*ny)+1]-2*ATP_cyto[0+0*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+(0+1)*nx+0*(nx*ny)]+ATP_cyto[0+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[0+0*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+0*nx+(0+1)*(nx*ny)]+ATP_cyto[0+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+0*nx+0*(nx*ny)])/tau_ATP;
  J_ATP_D[0+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP;
  J_ATP_D[0+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+0*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+0*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+0*(nx*ny)])/tau_ATP;
  J_ATP_D[0+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+0*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+0*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[0+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+(nz-1)*(nx*ny)])/tau_ATP;
  J_ATP_D[(nx-1)+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+0*(nx*ny)])/tau_ATP;
  J_ATP_D[(nx-1)+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+(nz-1)*(nx*ny)])/tau_ATP;
  J_ATP_D[(nx-1)+0*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+(0+1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+0*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+0*(nx*ny)])/tau_ATP;
  J_ATP_D[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
    (ATP_cyto[(nx-1)+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP;
 
  // x fixed (varies j, then k)
#pragma omp parallel for
  for (int j=1;j<ny-1;j++)
  {
    J_ATP_D[0+j*nx+0*(nx*ny)]=(ATP_cyto[0+j*nx+0*(nx*ny)+1]+ATP_cyto[0+j*nx+0*(nx*ny)+1]-2*ATP_cyto[0+j*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+(j+1)*nx+0*(nx*ny)]+ATP_cyto[0+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[0+j*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+j*nx+(0+1)*(nx*ny)]+ATP_cyto[0+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[0+j*nx+0*(nx*ny)])/tau_ATP;
    J_ATP_D[(nx-1)+j*nx+0*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+0*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+0*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+(j+1)*nx+0*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+j*nx+(0+1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+0*(nx*ny)])/tau_ATP;
    J_ATP_D[0+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[0+j*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[0+j*nx+(nz-1)*(nx*ny)+1]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[0+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[0+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+(nz-1)*(nx*ny)])/tau_ATP;
    J_ATP_D[(nx-1)+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+(nz-1)*(nx*ny)])/tau_ATP;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nz-1;k++)
    {
      J_ATP_D[0+j*nx+k*(nx*ny)]=(ATP_cyto[0+j*nx+k*(nx*ny)+1]+ATP_cyto[0+j*nx+k*(nx*ny)+1]-2*ATP_cyto[0+j*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[0+(j+1)*nx+k*(nx*ny)]+ATP_cyto[0+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[0+j*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[0+j*nx+(k+1)*(nx*ny)]+ATP_cyto[0+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+j*nx+k*(nx*ny)])/tau_ATP;
      J_ATP_D[(nx-1)+j*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+j*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+j*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[(nx-1)+(j+1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[(nx-1)+j*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+j*nx+k*(nx*ny)])/tau_ATP;
    }
  }
 
  // y fixed (varies i, then k, then j)
#pragma omp parallel for
  for (int i=1;i<(nx-1);i++)
  {
    J_ATP_D[i+0*nx+0*(nx*ny)]=(ATP_cyto[i+0*nx+0*(nx*ny)+1]+ATP_cyto[i+0*nx+0*(nx*ny)-1]-2*ATP_cyto[i+0*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+(0+1)*nx+0*(nx*ny)]+ATP_cyto[i+(0+1)*nx+0*(nx*ny)]-2*ATP_cyto[i+0*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+0*nx+(0+1)*(nx*ny)]+ATP_cyto[i+0*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+0*nx+0*(nx*ny)])/tau_ATP;
    J_ATP_D[i+(ny-1)*nx+0*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+0*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+0*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+((ny-1)-1)*nx+0*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+0*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+(ny-1)*nx+(0+1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+0*(nx*ny)])/tau_ATP;
    J_ATP_D[i+0*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+0*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+0*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+(0+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+(0+1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+0*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+0*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+(nz-1)*(nx*ny)])/tau_ATP;
    J_ATP_D[i+(ny-1)*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+((ny-1)-1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP+
      (ATP_cyto[i+(ny-1)*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+(nz-1)*(nx*ny)])/tau_ATP;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nz-1;k++)
    {
      J_ATP_D[i+0*nx+k*(nx*ny)]=(ATP_cyto[i+0*nx+k*(nx*ny)+1]+ATP_cyto[i+0*nx+k*(nx*ny)-1]-2*ATP_cyto[i+0*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+(0+1)*nx+k*(nx*ny)]+ATP_cyto[i+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[i+0*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+0*nx+(k+1)*(nx*ny)]+ATP_cyto[i+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+0*nx+k*(nx*ny)])/tau_ATP;
      J_ATP_D[i+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[i+(ny-1)*nx+k*(nx*ny)+1]+ATP_cyto[i+(ny-1)*nx+k*(nx*ny)-1]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[i+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[i+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+(ny-1)*nx+k*(nx*ny)])/tau_ATP;
    }
    // z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<ny-1;j++)
    {
      J_ATP_D[i+j*nx+0*(nx*ny)]=(ATP_cyto[i+j*nx+0*(nx*ny)+1]+ATP_cyto[i+j*nx+0*(nx*ny)-1]-2*ATP_cyto[i+j*nx+0*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+(j+1)*nx+0*(nx*ny)]+ATP_cyto[i+(j-1)*nx+0*(nx*ny)]-2*ATP_cyto[i+j*nx+0*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+j*nx+(0+1)*(nx*ny)]+ATP_cyto[i+j*nx+(0+1)*(nx*ny)]-2*ATP_cyto[i+j*nx+0*(nx*ny)])/tau_ATP;
      J_ATP_D[i+j*nx+(nz-1)*(nx*ny)]=(ATP_cyto[i+j*nx+(nz-1)*(nx*ny)+1]+ATP_cyto[i+j*nx+(nz-1)*(nx*ny)-1]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+(j+1)*nx+(nz-1)*(nx*ny)]+ATP_cyto[i+(j-1)*nx+(nz-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])/tau_ATP+
        (ATP_cyto[i+j*nx+((nz-1)-1)*(nx*ny)]+ATP_cyto[i+j*nx+((nz-1)-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+(nz-1)*(nx*ny)])/tau_ATP;
    }
  }
 
  // remaining k-fixed edges
#pragma ivdep
#pragma vector always
  for (int k=1;k<nz-1;k++)
  {
    J_ATP_D[0+0*nx+k*(nx*ny)]=(ATP_cyto[0+0*nx+k*(nx*ny)+1]+ATP_cyto[0+0*nx+k*(nx*ny)+1]-2*ATP_cyto[0+0*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+(0+1)*nx+k*(nx*ny)]+ATP_cyto[0+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[0+0*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+0*nx+(k+1)*(nx*ny)]+ATP_cyto[0+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+0*nx+k*(nx*ny)])/tau_ATP;
    J_ATP_D[0+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[0+(ny-1)*nx+k*(nx*ny)+1]+ATP_cyto[0+(ny-1)*nx+k*(nx*ny)+1]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[0+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[0+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[0+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[0+(ny-1)*nx+k*(nx*ny)])/tau_ATP;
    J_ATP_D[(nx-1)+0*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+0*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+0*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+(0+1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+(0+1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+0*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+0*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+0*nx+k*(nx*ny)])/tau_ATP;
    J_ATP_D[(nx-1)+(ny-1)*nx+k*(nx*ny)]=(ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)-1]+ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)-1]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+((ny-1)-1)*nx+k*(nx*ny)]+ATP_cyto[(nx-1)+((ny-1)-1)*nx+k*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])/tau_ATP+
      (ATP_cyto[(nx-1)+(ny-1)*nx+(k+1)*(nx*ny)]+ATP_cyto[(nx-1)+(ny-1)*nx+(k-1)*(nx*ny)]-2*ATP_cyto[(nx-1)+(ny-1)*nx+k*(nx*ny)])/tau_ATP;
  }
  
  // interior — covers the vast majority of the grid
#pragma omp parallel for
  for (int k=1;k<nz-1;k++)
  {
    for (int j=1;j<ny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nx-1;i++)
      {
        J_ATP_D[i+j*nx+k*(nx*ny)]=(ATP_cyto[i+j*nx+k*(nx*ny)+1]+ATP_cyto[i+j*nx+k*(nx*ny)-1]-2*ATP_cyto[i+j*nx+k*(nx*ny)])/tau_ATP+
          (ATP_cyto[i+(j+1)*nx+k*(nx*ny)]+ATP_cyto[i+(j-1)*nx+k*(nx*ny)]-2*ATP_cyto[i+j*nx+k*(nx*ny)])/tau_ATP+
          (ATP_cyto[i+j*nx+(k+1)*(nx*ny)]+ATP_cyto[i+j*nx+(k-1)*(nx*ny)]-2*ATP_cyto[i+j*nx+k*(nx*ny)])/tau_ATP;
      }
    }
  }
}   
